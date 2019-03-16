#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/spirit/include/support_line_pos_iterator.hpp>
#include <programOptions.h>
#include <unordered_map>
#include <fstream>

namespace gtf{
   struct TransFrag{
      std::string feature;
      std::string score;
      std::string frame;
      int start;
      int end;
   };

   struct TransFragFull:TransFrag{
      std::string seqname;
      std::string source;
      std::string strand;
      std::unordered_map<std::string,std::string> attributes;
   };

   struct Transcript{
      std::string seqname;
      int chrVal;
      std::string source;
      std::string strand;
      std::string transcript_type;
      std::string transcript_name;
      std::string gene_id;
      std::string gene_type;
      std::string gene_name;
      std::vector<TransFrag> fragments;
   };

   typedef std::unordered_map<std::string,Transcript> TranscriptMap;

   template<typename Iter>
   struct transfrag_parser : boost::spirit::qi::grammar<Iter, TranscriptMap()>{
      static bool addToMap(TranscriptMap &fragMap, TransFragFull &frag, const std::string &suffix){
         const std::string transKey=frag.attributes["transcript_id"]+suffix; //custom suffix, e.g. used to distinguish between transcripts from different files
         auto mapIter=fragMap.find(transKey);
         if(mapIter==fragMap.end()){
            if(NUMCHR>0){
               if(frag.seqname.compare(0,3,"chr")){
                  warn("Invalid chromosome ("+frag.seqname+")");
                  return 0;
               }
               const std::string chr = frag.seqname.substr(3);
               const int bigVal=std::numeric_limits<int>::max(); //X, Y and M chromosomes should be ordered last
               if(chr=="X") fragMap[transKey].chrVal=bigVal-2;
               else if(chr=="Y") fragMap[transKey].chrVal=bigVal-1;
               else if(chr== "M") fragMap[transKey].chrVal=bigVal;
               else{
                  fragMap[transKey].chrVal=atoi(chr.c_str());
                  if(fragMap[transKey].chrVal<1||fragMap[transKey].chrVal>NUMCHR){
                     warn("Invalid chromosome: "+frag.seqname);
                     return 0;
                  }
               }
            }
            fragMap[transKey].seqname=frag.seqname;
            fragMap[transKey].source=frag.source;
            fragMap[transKey].strand=frag.strand;
            fragMap[transKey].transcript_type=frag.attributes["transcript_type"];
            fragMap[transKey].transcript_name=frag.attributes["transcript_name"];
            fragMap[transKey].gene_id=frag.attributes["gene_id"]+suffix;
            fragMap[transKey].gene_type=frag.attributes["gene_type"];
            fragMap[transKey].gene_name=frag.attributes["gene_name"];
            frag.attributes.clear();
         }
         else{
            if(frag.strand!=mapIter->second.strand){
               fragMap["strand\n"]; return 1;} //extremely hacky way of indicating what failed
            if(frag.seqname!=mapIter->second.seqname){ //while preserving line number in exception
               fragMap["seqname\n"]; return 1;}
            if(frag.strand=="."){
               fragMap["multidot\n"]; return 1;}
         }
         std::vector<TransFrag> &trans=fragMap[transKey].fragments;
         auto transIter=trans.rbegin(); //add fragments in order of genomic start position
         for(; transIter!=trans.rend(); ++transIter) if(transIter->start<frag.start) break;
         trans.insert(transIter.base(),frag);
         return 0;
      }

      //helper function to print parser warning message
      static void warn(const std::string &msg){
         BOOST_LOG_TRIVIAL(warning)<<"Warning: "<<msg;}

      //helper function to store attributes from parsing
      static void setAttribute(TransFragFull &frag,std::pair<std::string,std::string> attribute){
         frag.attributes[attribute.first]=attribute.second;}

      transfrag_parser(const std::string &suffix) : transfrag_parser::base_type(start){
         using boost::spirit::qi::int_;
         using boost::spirit::qi::omit;
         using boost::spirit::ascii::char_;
         using boost::phoenix::bind;
         using boost::spirit::_val;
         using boost::spirit::_1;
         using boost::spirit::eol;
         using boost::spirit::eoi;
         using boost::spirit::ascii::blank;
         using boost::spirit::qi::graph;
         using boost::spirit::qi::_pass;
         using boost::spirit::lit;

         string_field %= +graph;
         pair %= key > value > -blank;
         key %= *(graph-' '-'\"');
         value %= ((lit(" \"")|lit("\"")) > *(char_-'\"') > "\";") | (' '>*(graph-';') > ';');
         comment %= ('#'|(+blank > '#')) > *(char_-eol-eoi);

         frag %=
            string_field[bind(&TransFragFull::seqname,_val)=_1] > '\t'
            > string_field[bind(&TransFragFull::source,_val)=_1] > '\t'
            > string_field[bind(&TransFragFull::feature,_val)=_1] > '\t'
            > int_[bind(&TransFragFull::start,_val)=_1] > '\t'
            > int_[bind(&TransFragFull::end,_val)=_1] > '\t'
            > string_field[bind(&TransFragFull::score,_val)=_1] > '\t'
            > string_field[bind(&TransFragFull::strand,_val)=_1] > '\t'
            > string_field[bind(&TransFragFull::frame,_val)=_1] > '\t'
            > omit[*(pair[bind(&setAttribute,_val,_1)]-blank-'#'-eol-eoi)] > -comment
            ;
            start %= *(omit[(comment|frag[_pass=!bind(&addToMap,_val,_1,suffix)]) > (eol|eoi)])
               > omit[eoi|(*eol>eoi)[bind(&warn,"There are blank lines at the end of the GTF")]];
      }

      boost::spirit::qi::rule<Iter, std::string()> string_field, key, value, comment, blanks;
      boost::spirit::qi::rule<Iter, std::pair<std::string,std::string>()> pair;
      boost::spirit::qi::rule<Iter, TransFragFull()> frag;
      boost::spirit::qi::rule<Iter, TranscriptMap()> start;
   };

   class gtf{
   public:
      gtf(){};
      gtf(const std::string &fName, const std::string &suffix=""){
         readFromFile(fName,suffix);
      }

      void readFromFile(const std::string &fName, const std::string &suffix=""){
         boost::iostreams::mapped_file_source mmap(fName.c_str()); //memory mapped for speed
         typedef boost::spirit::line_pos_iterator<boost::iostreams::mapped_file::const_iterator> line_pos_iterator_type;
         line_pos_iterator_type position_begin(mmap.data());
         line_pos_iterator_type position_end(mmap.end());
         transfrag_parser<line_pos_iterator_type> g(suffix);
         try{
            parse(position_begin, position_end, g, transcripts);
         } catch(const boost::spirit::qi::expectation_failure<line_pos_iterator_type>& e){
            std::string message=e.what();
            if(transcripts.count("strand\n")) message="mismatched strand";
            else if(transcripts.count("seqname\n")) message="mismatched sequence name";
            else if(transcripts.count("multidot\n")) message="multiple fragments with '.' strand";
            throw std::runtime_error("Error: problem reading file at line "+std::to_string(get_line(e.first))+" ("+message+")\n");
         }
      }

      void writeToFile(const std::string &fName, bool append=false){
         std::ofstream outFile;
         if(append)
            outFile.open(fName.c_str(), std::ofstream::out | std::ofstream::app);
         else
            outFile.open(fName.c_str());

         if(NUMCHR>0)
            writeToFile<int>(outFile);
         else
            writeToFile<std::string>(outFile);
         outFile.close();
      }
      TranscriptMap transcripts;

   private:
      void getChr(int &chr, const Transcript &trans){chr=trans.chrVal;}
      void getChr(std::string &chr, const Transcript &trans){chr=trans.seqname;}

      template<typename T>
      void writeToFile(std::ofstream &outFile){
         std::multimap<std::pair<T,int>,std::string> outMap;
         for (auto &trans : transcripts){
            int exonIdx=0;
            for(auto &frag : trans.second.fragments){
               std::ostringstream line;
               line << trans.second.seqname << '\t' << trans.second.source;
               line << '\t' << frag.feature << '\t' << frag.start;
               line << '\t' << frag.end << '\t' << frag.score;
               line << '\t' << trans.second.strand << '\t' << frag.frame << '\t';
               line << "exon_number \"" << ++exonIdx << "\"; ";
               line << "gene_id \"" << trans.second.gene_id << "\"; ";
               line << "gene_type \"" << trans.second.gene_type << "\"; ";
               line << "gene_name \"" << trans.second.gene_name << "\"; ";
               line << "transcript_id \"" << trans.first << "\"; ";
               line << "transcript_type \"" << trans.second.transcript_type << "\"; ";
               line << "transcript_name \"" << trans.second.transcript_name << "\";\n";
               T chr; getChr(chr,trans.second);
               outMap.insert({{chr,frag.start},line.str()});
            }
         }
         for(auto &line : outMap)
            outFile << line.second;
      }
   };
}

BOOST_FUSION_ADAPT_STRUCT(gtf::TransFragFull,seqname,source,feature,start,end,score,strand,frame)
