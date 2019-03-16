#include <boost/log/trivial.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>
#define STR(X) static_cast<std::ostringstream&>(std::ostringstream()<<X).str()

namespace gtf{
   namespace po = boost::program_options;
   int NUMCHR=22; //used to test for invalid chromosome number

   struct base_option{ //to be inherited from
      base_option(const std::string &name, const std::string &desc):name(name),desc(desc){};
      virtual void add(po::options_description &options){}
      virtual void process(po::variable_value poVal){}
      virtual std::string text(){}
      virtual bool isRequired(){return false;}
      const std::string name;
      const std::string desc;
   };

   struct flag:base_option{ //option without value (i.e. flag)
      flag(const std::string &name, const std::string &desc, std::function<void()> const &action) :
         base_option(name,desc),action(action){};
      void add(po::options_description &options){
         options.add_options()(name.c_str(),desc.c_str());}
      std::string text(){if(set)return desc+" (--"+name+")"; else return "";}
      void process(po::variable_value poVal){action();set=true;}
      bool set=false;
      std::function<void()> const action;
   };

   template <typename T>
   struct option:base_option{ //option with value (input argument)
      option(const std::string name, const std::string &desc, const bool required, T *val) :
         base_option(name,desc),required(required),val(val){};
      void add(po::options_description &options){
         options.add_options()(name.c_str(),po::value<T>(),
                               STR(desc<<" (default: "<<*val<<+")").c_str());
      }
      void process(po::variable_value poVal){*val=poVal.as<T>();}
      std::string text(){return STR(desc<<" (--"<<name<<"): "<<*val);}
      bool isRequired(){return required;}
      const bool required;
      T *val;
   };

   class programOptions{
   public:
      programOptions(const std::string &programTitle="",
                     const bool fileOutput=false,
                     const int outputFiles=-1) : programTitle(programTitle),
                                                 fileOutput(fileOutput),
                                                 inputFiles(inputFiles){};

      ~programOptions(){
         const std::time_t endTime = std::time(nullptr);
         BOOST_LOG_TRIVIAL(info)<<"End time: "<<std::asctime(std::localtime(&endTime));
      }

      template <typename T>
      void addOption(T &val, const std::string &name="", const std::string &description="", const bool required=false, const T &defaultValue=0){
         options.push_back(std::unique_ptr<option<T>>(new option<T>(name,description,required,&(val=defaultValue))));
      }

      void addFlag(const std::string &name="", const std::string &description="", std::function<void()> const &action=std::function<void()>()){
         options.push_back(std::unique_ptr<flag>(new flag(name,description,action)));
      }

      void parse(int argc, const char* argv[]){
         const std::time_t startTime = std::time(nullptr);
         boost::log::add_console_log(std::cout);
         po::options_description desc("Command line options");

         addFlag("help","Produce help message",[&desc](){BOOST_LOG_TRIVIAL(info)<<desc;exit(0);});
         addOption(NUMCHR,"numChromosomes","Number (N) of numeric autosomes",false,22);
         addFlag("allowOtherSequenceNames","Do not restrict to chr[1-N,X,Y,M]",[](){NUMCHR=0;});
         if(fileOutput)
            addOption(outPrefix,"outPrefix","Set output prefix (default: \"./out\")",false,std::string("./out"));
         for(auto &opts : options)
            opts->add(desc);

         po::options_description desc_hidden("Hidden options"); //do not want --files to appear in help
         desc_hidden.add_options()("files", po::value(&file_names), "list of files");
         po::options_description cmdline;
         cmdline.add(desc).add(desc_hidden);
         po::positional_options_description pos; //do not need to specify --files
         pos.add("files", -1);
         
         po::variables_map vm;
         try {
            if(inputFiles>=0)
               po::store(po::command_line_parser(argc,argv).options(cmdline).positional(pos).run(),vm);
            else
               po::store(po::command_line_parser(argc,argv).options(desc).run(),vm);
         }
         catch(const po::error& e) {
            throw std::runtime_error(STR("Error: Could not parse command line arguments properly:\n"<<e.what()<<"\n"<<desc));
         }
         po::notify(vm);

         for(auto &opts : options)
            if(vm.count(opts->name)) opts->process(vm[opts->name]);
            else if(opts->isRequired()) throw std::runtime_error("Error: --"+opts->name+" required, but not provided");

         if(fileOutput){
            if(vm.count("outPrefix")) outPrefix="./"+outPrefix;
            boost::log::add_file_log(outPrefix+".log"); //create log of parameters used
         }

         if(file_names.size()<inputFiles){
            throw std::runtime_error(STR("Error: gtfmerge requires at least 2 GTF files to merge\n"<<desc));
         }

         BOOST_LOG_TRIVIAL(info)<<programTitle<<"\n";
         BOOST_LOG_TRIVIAL(info)<<"Options in effect:";
         for(auto &opts : options){
            std::string text=opts->text();
            if(!text.empty()) BOOST_LOG_TRIVIAL(info)<<"\t"<<text;
         }

         if(inputFiles>=0)
            for (int fIdx=0; fIdx<file_names.size(); ++fIdx)
               BOOST_LOG_TRIVIAL(info)<<"\tGTF "<<fIdx<<": "<<file_names[fIdx];

         BOOST_LOG_TRIVIAL(info)<<"\nWorking directory: "<<boost::filesystem::current_path();
         BOOST_LOG_TRIVIAL(info)<<"Start time: "<<std::asctime(std::localtime(&startTime));
      }

      std::string outPrefix="./out";
      std::vector<std::string> file_names;
   private:
      const std::string programTitle;
      const bool fileOutput;
      const int inputFiles;
      std::vector<std::unique_ptr<base_option>> options;
   };
}
