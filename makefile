gtfmerge: gtfMerge.cpp gtf.h
  g++ -std=c++11 -O2 -o gtfmerge -I. gtfMerge.cpp -lboost_iostreams -lboost_program_options -lboost_system -lboost_filesystem

clean:
  rm gtfmerge
