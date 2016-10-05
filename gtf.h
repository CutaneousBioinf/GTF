include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/spirit/include/support_line_pos_iterator.hpp>
#include <unordered_map>
#include <fstream>

namespace gtf{

   struct transfrag{
      std::string seqname;
      std::string source;
      std::string feature;
      int start;
      int end;
      std::string score;
      std::string strand;
      std::string frame;
      std::string transcript_id;
      std::string transcript_name;
      std::string gene_id;
      std::string gene_type;
      std::string gene_name;
      std::unordered_map<std::string,std::string> attributes;
   };

   typedef std::unordered_map<std::string,std::vector<transfrag>> TranscriptMap;

   template<typename Iter,char suffix>
   struct transfrag_parser : boost::spirit::qi::grammar<Iter, TranscriptMap()>{
      static void addToMap(TranscriptMap &fragMap,transfrag &frag){
        frag.transcript_id=frag.attributes["transcript_id"];
        frag.transcript_type=frag.attributes["transcript_type"];
        frag.transcript_name=frag.attributes["transcript_name"];
        frag.gene_id=frag.attributes["gene_id"];
        frag.gene_type=frag.attributes["gene_type"];
        frag.gene_name=frag.attributes["gene_name"];
        frag.attributes.clear();
        std::vector<transfrag> &trans=fragMap[frag.transcript_id+suffix];
        auto transIter=trans.rbegin();
        for(; transIter!=trans.rend(); ++transIter) if(transIter->start<frag.start) break;
        trans.insert(transIter.base(),frag);
      }
      static void setAttribute(transfrag &frag,std::pair<std::string,std::string> attribute){
         frag.attributes[attribute.first]=attribute.second;}

      transfrag_parser() : transfrag_parser::base_type(start){
         using boost::spirit::qi::int_;
         using boost::spirit::qi::lexeme;
         using boost::spirit::qi::omit;
         using boost::spirit::qi::lit;
         using boost::spirit::ascii::char_;
         using boost::phoenix::bind;
         using boost::spirit::_val;
         using boost::spirit::_1;
         using boost::spirit::eol;
         using boost::spirit::ascii::blank;

         string_field %= +(char_ - '\t');
         att_arg %= " \"">*(char_ - '\"')>"\";";
         pair %= key > value > -lit(" ");
         key %= *(char_-" ")>" ";
         value %= ("\"">*(char_-"\""-eol)>"\";") | (*(char_-";"-eol)>";");

         frag %=
            string_field[bind(&transfrag::seqname,_val)=_1] > '\t'
            > string_field[bind(&transfrag::source,_val)=_1] > '\t'
            > string_field[bind(&transfrag::feature,_val)=_1] > '\t'
            > int_[bind(&transfrag::start,_val)=_1] > '\t'
            > int_[bind(&transfrag::end,_val)=_1] > '\t'
            > string_field[bind(&transfrag::score,_val)=_1] > '\t'
            > string_field[bind(&transfrag::strand,_val)=_1] > '\t'
            > string_field[bind(&transfrag::frame,_val)=_1] > '\t'
            > omit[*(pair[bind(&setAttribute,_val,_1)] - '#' - eol - '\t' - "  ")] > *blank > -('#'>*(char_-eol)) > eol
            ;
            start %= *(omit[*('#'>*(char_-eol)>eol)>>frag[bind(&addToMap,_val,_1)]]);
      }

      boost::spirit::qi::rule<Iter, std::string()> string_field, att_arg, key, value;
      boost::spirit::qi::rule<Iter, std::pair<std::string,std::string>()> pair;
      boost::spirit::qi::rule<Iter, transfrag()> frag;
      boost::spirit::qi::rule<Iter, TranscriptMap()> start;
   };

   template <char suffix>
   class gtf{
   public:
      gtf(const std::string &fName){readFromFile(fName);}
      gtf(){};
      void readFromFile(const std::string &fName){
         boost::iostreams::mapped_file_source mmap(fName.c_str());
         typedef boost::spirit::line_pos_iterator<boost::iostreams::mapped_file::const_iterator> line_pos_iterator_type;
         line_pos_iterator_type position_begin(mmap.data());
         line_pos_iterator_type position_end(mmap.end());
         transfrag_parser<line_pos_iterator_type,suffix> g;
         try{
            parse(position_begin, position_end, g, transcripts);
         } catch(const boost::spirit::qi::expectation_failure<line_pos_iterator_type>& e){
            std::cout<<"Error: problem reading file at line "<<get_line(e.first)<<std::endl;
            exit(EXIT_FAILURE);
         }
      }

   struct keyComparator
   {
      keyComparator(TranscriptMap &transcripts) : transcripts(transcripts) {}
      inline bool operator() (const std::string& key1, const std::string& key2)
      {
         transfrag &tf1=transcripts[key1][0], &tf2=transcripts[key2][0];
         if(tf1.seqname==tf2.seqname)
            return (tf1.start < tf2.start);
         else
            return (tf1.seqname < tf2.seqname);
      }
      TranscriptMap &transcripts;
   };

      void writeToFile(const std::string &fName, bool append=false){
        std::vector<std::string> keys;
        keys.reserve (transcripts.size());
        for (auto& it : transcripts) {
           keys.push_back(it.first);
        }
        std::sort(keys.begin(),keys.end(),keyComparator(transcripts));

        std::ofstream outFile;
        if(append)
            outFile.open(fName.c_str(), std::ofstream::out | std::ofstream::app);
         else
            outFile.open(fName.c_str());
         for(auto keyIter=keys.begin(); keyIter!=keys.end(); ++keyIter){
            char transcript_suffix = keyIter->at(keyIter->size()-1);
            int exonIdx=0;
            for(auto fragIter=transcripts[*keyIter].begin(); fragIter!=transcripts[*keyIter].end(); ++fragIter){
               outFile << fragIter->seqname << '\t' << fragIter->source << '\t' << fragIter->feature << '\t';
               outFile << fragIter->start << '\t' << fragIter->end;
               outFile << '\t' << fragIter->score << '\t' << fragIter->strand << '\t' << fragIter->frame << '\t';
               outFile << "exon_number \"" << ++exonIdx << "\"; ";
               outFile << "gene_id \"" << fragIter->gene_id;
               outFile << "_" << transcript_suffix << "\"; ";
               outFile << "gene_type \"" << fragIter->gene_type << "\"; ";
               outFile << "gene_name \"" << fragIter->gene_name << "\"; ";
               outFile << "transcript_id \"" << fragIter->transcript_id;
               outFile << "_" << transcript_suffix << "\"; ";
               outFile << "transcript_type \"" << fragIter->transcript_type << "\"; ";
               outFile << "transcript_name \"" << fragIter->transcript_name << "\";\n";
            }
         }
         outFile.close();
      }
      TranscriptMap transcripts;
   };
}

BOOST_FUSION_ADAPT_STRUCT(gtf::transfrag,)
