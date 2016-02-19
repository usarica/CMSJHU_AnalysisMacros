#include "../interface/Reader.h"
#include "../interface/convertLHE.h"
#include "../interface/convertPythia.h"

int main(int argc, char** argv){
  OptionParser options(argc, argv);
  options.printOptionSummary();
  if (options.analysisLevel()==-1) Reader analyzer(&options);
  else if (options.analysisLevel()==0) convertLHE analyzer(&options);
  else if (options.analysisLevel()==1 || options.analysisLevel()==2) convertPythia analyzer(&options);
}

