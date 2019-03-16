gtfmerge: gtfMerge.cpp gtf.h
        g++ -std=c++11 -O2 -static -o gtfmerge -I. gtfMerge.cpp -lboost_log -lboost_system -lboost_thread -lpthread -lboost_iostreams -lboost_filesystem -lboost_program_options

clean:
        rm gtfmerge
