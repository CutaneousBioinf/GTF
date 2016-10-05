#include <gtf.h>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

namespace gtf{

   bool isOverlapping(std::vector<transfrag> &tf1, std::vector<transfrag> &tf2, int minBasePairs=50, double minPercentage=20.0){
      int count=0, fragmentSize=0; auto tf2Iter=tf2.begin();
      for(auto tf1Iter=tf1.begin(); tf1Iter!=tf1.end(); ++tf1Iter){
         fragmentSize+=tf1Iter->end-tf1Iter->start;
         while(tf2Iter!=tf2.end()){
            if(tf2Iter->end > tf1Iter->start){break;}
            tf2Iter++;
         }
         while(tf2Iter!=tf2.end()){
            int c=std::min(tf2Iter->end,tf1Iter->end)-std::max(tf2Iter->start,tf1Iter->start)+1;
            if(c > 0){count+=c;} else{break;}
            tf2Iter++;
         }
         if(tf2Iter==tf2.end()){break;} else{tf2Iter--;}
      }
      if(count>=minBasePairs || count>=minPercentage*fragmentSize){return true;} else {return false;}
   }

   void merge(TranscriptMap &ref, TranscriptMap &target, TranscriptMap &add, int minBasePairs, double minPercentage){
      for(auto targetIter=target.begin(); targetIter!=target.end();){
         bool overlap=false;
         for(auto refIter=ref.begin(); refIter!=ref.end(); ++refIter){
            if(targetIter->second.begin()->seqname == refIter->second.begin()->seqname)
               if(targetIter->second.begin()->strand==refIter->second.begin()->strand || targetIter->second.begin()->strand=="." || refIter->second.begin()->strand==".")
                  if(isOverlapping(targetIter->second,refIter->second,minBasePairs,minPercentage)){overlap=true; break;}
         }
         if(!overlap){
            add[targetIter->first]=targetIter->second;
            targetIter=target.erase(targetIter);
         }
         else
            targetIter++;
      }
      for (auto addIter=add.begin(); addIter!=add.end(); ++addIter){
         ref[addIter->first]=addIter->second;
      }
   }
}

int main(int argc, const char* argv[])
{
   namespace po = boost::program_options;
   po::options_description desc("Command line options");
   desc.add_options()
      ("help","produce help message")
      ("outDir",po::value<std::string>(),"set output directory (default: \"./\")")
      ("outPre",po::value<std::string>(),"set output prefix (default: \"\")")
      ("minBP",po::value<int>(),"set minimum number basepairs for overlap (default: 50)")
      ("minPct",po::value<double>(),"set minimum percentage for overlap (default: 20.0)")
      ("refGTF",po::value<std::string>(),"set reference file (can be specified without the flag)")
      ("targetGTF",po::value<std::string>(),"set target file (can be specified without the flag)")
   ;

   po::positional_options_description p;
   p.add("refGTF",1).add("targetGTF",1);

   po::variables_map vm;
   po::store(po::command_line_parser(argc,argv).options(desc).positional(p).run(),vm);
   po::notify(vm);

   if(vm.count("help")){
      std::cout<<desc<<std::endl;
      exit(0);
   }

   std::string outDir=".",outPrefix="",ref="",target="";
   int minBasePairs=50; double minPercentage=20.0;
   if(vm.count("outDir")) outDir=vm["outDir"].as<std::string>();
   if(vm.count("outPre")) outPrefix=vm["outPre"].as<std::string>();
   if(vm.count("minBP")) minBasePairs=vm["minBP"].as<int>();
   if(vm.count("minPct")) minPercentage=vm["minPct"].as<int>();
   if(vm.count("refGTF")) ref=vm["refGTF"].as<std::string>();
   else{std::cerr<<"Error: reference file must be specified"<<std::endl; exit(EXIT_FAILURE);}
   if(vm.count("targetGTF")) target=vm["targetGTF"].as<std::string>();
   else{std::cerr<<"Error: target file must be specified"<<std::endl; exit(EXIT_FAILURE);}

   boost::filesystem::path outPath(outDir);
   if(!boost::filesystem::is_directory(outDir)){
      std::cout<<"Warning:: creating output directory ("<<outDir<<") as it does not currently exist"<<std::endl;
      boost::filesystem::create_directory(outDir);
   }

   std::cout<<"Reading target..."<<std::endl;
   gtf::gtf<'1'> g2(target);
   std::cout<<"Reading reference..."<<std::endl;
   gtf::gtf<'0'> g1(ref); gtf::gtf<'_'> g3;
   std::cout<<"Merging target and reference..."<<std::endl;

   gtf::merge(g1.transcripts,g2.transcripts,g3.transcripts,minBasePairs,minPercentage);
   std::cout<<"Writing output..."<<std::endl;
   g1.writeToFile(outDir+"/"+outPrefix+"combined.gtf");
   g2.writeToFile(outDir+"/"+outPrefix+"overlap.gtf");
   g3.writeToFile(outDir+"/"+outPrefix+"added.gtf");
}            
