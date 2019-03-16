#include <gtf.h>

namespace gtf{

   bool isOverlapping(Transcript &tf1, Transcript &tf2, int minBasePairs=50, double minPercentage=20.0){
      int count=0, fragmentSize=0; auto tf2Iter=tf2.fragments.begin();
      for(auto tf1Iter=tf1.fragments.begin(); tf1Iter!=tf1.fragments.end(); ++tf1Iter){
         fragmentSize+=tf1Iter->end-tf1Iter->start+1;
         while(tf2Iter!=tf2.fragments.end()){ //find nearest tf2 to tf1
            if(tf2Iter->end > tf1Iter->start){break;}
            tf2Iter++;
         }
         while(tf2Iter!=tf2.fragments.end()){
            int c=std::min(tf2Iter->end,tf1Iter->end)-std::max(tf2Iter->start,tf1Iter->start)+1;
            if(c > 0){count+=c;} else{break;} //positive values of c indicate some overlap
            tf2Iter++;
         }
         if(tf2Iter==tf2.fragments.end()) break; //no more tf2s to overlap with tf1s
      }
      if(count>=minBasePairs || count>=minPercentage*fragmentSize){return true;} else {return false;}
   }

   void merge(TranscriptMap &ref, TranscriptMap &target, TranscriptMap &add, int minBasePairs, double minPercentage){
      for(auto targetIter=target.begin(); targetIter!=target.end();){
         bool overlap=false; //loop through all transcripts each time, as TranscriptMap is unordered
         for(auto refIter=ref.begin(); refIter!=ref.end(); ++refIter){
            if(targetIter->second.seqname == refIter->second.seqname)
               if(targetIter->second.strand==refIter->second.strand || targetIter->second.strand=="." || refIter->second.strand==".") //transcripts containing only one fragment may be in either order
                  if(isOverlapping(targetIter->second,refIter->second,minBasePairs,minPercentage)){overlap=true; break;}
         }
         if(!overlap){
            add[targetIter->first]=targetIter->second;
            targetIter=target.erase(targetIter); //points to next element
         }
         else
            targetIter++;
      }
      for (auto addIter=add.begin(); addIter!=add.end(); ++addIter){
         ref[addIter->first]=addIter->second; //important to add at the end, to avoid interfering with loop
      }
   }
}

int main(int argc, const char* argv[])
{
   try{
      int minBasePairs;
      double minPercentage;
      gtf::programOptions opt("gtfmerge v0.1 (9 March 2019)",true,2);
      opt.addOption(minBasePairs,"minBP","Minimum number of basepairs for overlap",true,50);
      opt.addOption(minPercentage,"minPercent","Minimum percentage for overlap)",true,20.0);
      opt.parse(argc,argv);

      if(opt.file_names.size()<2){
         BOOST_LOG_TRIVIAL(error) << "Error: gtfmerge requires at least 2 GTF files to merge";
         return(1);
      }

      BOOST_LOG_TRIVIAL(info)<<"Reading GTF 0...";
      gtf::gtf g1(opt.file_names[0],"_0");  //first reference file
      for(int addIdx=1; addIdx<opt.file_names.size(); ++addIdx){
         BOOST_LOG_TRIVIAL(info)<<"Reading GTF "<<addIdx<<"...";
         gtf::gtf g2(opt.file_names[addIdx],"_"+std::to_string(addIdx));
         gtf::gtf g3; //temporary GTF to store added transcripts
         BOOST_LOG_TRIVIAL(info)<<"Merging GTF "<<addIdx-1<<" and "<<addIdx<<"...";
         gtf::merge(g1.transcripts, //reference GTF from first round becomes reference for next
                    g2.transcripts, //target GTF from first round becomes temporary GTF for overlap
                    g3.transcripts,
                    minBasePairs,minPercentage);
         BOOST_LOG_TRIVIAL(info)<<"Writing overlapping and added transcripts from GTF "<<addIdx<<"...\n";
         g2.writeToFile(opt.outPrefix+"_overlap"+std::to_string(addIdx)+".gtf");
         g3.writeToFile(opt.outPrefix+"_added"+std::to_string(addIdx)+".gtf");
      }
      BOOST_LOG_TRIVIAL(info)<<"Writing output...\n";
      g1.writeToFile(opt.outPrefix+"_combined.gtf");
   }
   catch(const std::exception &e){
      BOOST_LOG_TRIVIAL(error)<<e.what();
   }
}
