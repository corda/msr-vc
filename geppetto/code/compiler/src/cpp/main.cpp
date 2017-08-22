/*
 *  Run with -h for all of the options.  Also checks config.txt for option specifications.
 */
#pragma warning( disable : 4005 )	// Ignore warning about intsafe.h macro redefinitions
#include "Types.h"
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include "timer.h"

using namespace std;
using namespace boost::program_options;

TimerList timers;
extern void microbenchmarks(variables_map& config);
extern void functionality_tests(variables_map& config);

// Does the real work of main(), so memory-leak detection will be simpler
// (returning from here should clear everything off the stack)
void dispatch(int argc, char **argv) {
	// Parse the command-line options
	options_description genericOpt("Generic options");
	genericOpt.add_options()
			("help,h", "Produce the help message")
			("config,c", value<string>(), "Specify a config file.  Defaults to config.txt.")
			("raw", "Dump output in a raw, machine-readable format, instead of human readable")
			;

	options_description memOpt("Memory options");
	memOpt.add_options()				
			("mem", value<int>(), "REQUIRED: Maximum amount of memory to use for crypto tables (in GB)")
			("fill", "Precomputation is 'free', so fill all of memory with crypto tables")			
			;

	options_description highOpt("High-level options (must choose one)");
	highOpt.add_options()
			("qap", "Build a QAP")
			("micro", value<string>(), "Run microbenchmarks.  Options include: 'all', 'field', 'encoding', 'poly-all', 'poly-fast', 'poly-slow', 'mont-test'")
			("test", "Run basic tests of functionality.")
			;

	options_description encOpt("Additional options");
	encOpt.add_options()			
			("encoding,e", value<string>(), "REQUIRED. Which encoding to use.  Options include: 'enigma-bn', 'arith-bn', 'debug', 'arith-cp', 'arith-qap', 'arith-cp-dbg', 'arith-bn-dbg', 'arith-cp3'")			
			("trials,t", value<int>(), "REQUIRED for 'micro' and 'test'. Number of iterations to use.")			
			;

	// Group all of the options together
	options_description all("Allowed options");
	all.add(genericOpt).add(highOpt).add(memOpt).add(encOpt);

	variables_map config;
	store(parse_command_line(argc, argv, all), config);
	notify(config);  

	if (config.count("help")) {
		cout << all << endl;
		return;
	}

	string config_filename = "config.txt";		// Default value
	if (config.count("config")) {
		config_filename = config["config"].as<string>();
	} 

	ifstream ifs(config_filename.c_str());
  if (!ifs) {
      cout << "Warning: Could not open config file: " 
			     << config_filename 
					 << ".  Continuing based on command line only.\n";
  } else {
      store(parse_config_file(ifs, all), config);
      notify(config);
			ifs.close();
  }

	if (!config.count("mem")) {
		cout << "Error: Must specify how much memory to use for crypto tables\n";
		return;
	}

	if (!config.count("encoding")) {
		cout << "Error: Must specify which encoding to use.\n";
		return;
	}

	if ((config.count("micro") || config.count("test")) && !config.count("trials")) {
		cout << "Error: Need to tell me how many trials to use.\n";
		return;
	}

// Display info about the current configuration
#ifdef DEBUG_ENCODING
	printf("Using the identity encoding!\n");
#endif 

#ifdef DEBUG_RANDOM
	printf("Warning!  Random elements aren't that random!\n");
  srand(1402345247);
#else
	unsigned int seed = (unsigned int)time(NULL);	// Save this for debugging purposes
	srand(seed);
#endif 

#ifdef _DEBUG
	printf("\nRunning in DEBUG mode -- performance measurements should not be trusted!\n\n");
#endif

#ifdef USE_OLD_PROTO
	printf("\nUsing the OLD crypto protocol!\n");
#else !USE_OLD_PROTO
	printf("\nUsing the NEW crypto protocol!\n");
#endif

//#ifndef _DEBUG	// If we're debugging, let VS catch the exceptions, so we can see where they came from
//	try {
//#endif
		if (config.count("qap")) {
			//runQapTests(config);
		} else if (config.count("micro")) {
			microbenchmarks(config);
		} else if (config.count("test")) {
			functionality_tests(config);
		} else {
			printf("\nWarning: No high-level option recognized, so I'm not doing anything.\n");
		}
//#ifndef _DEBUG
//	} catch (exception e) {
//		cout << "Caught an exception: " << e.what() << endl;
//	} catch (...) {
//		cout << "Caught an unknown exception\n";
//	}
//#endif
	
	timers.cleanup();

}

#ifndef APPSDBG
int main(int argc, char **argv) {

#ifdef MEM_LEAK
  //_CrtSetBreakAlloc(1941);
#endif

	dispatch(argc, argv);

	#ifdef MEM_LEAK
	//_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
	//_CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);  
	_CrtDumpMemoryLeaks();  // Check for memory leaks
	#endif

	printf("\nDone!\n");
#ifdef _DEBUG
	getc(stdin); // Make the VS window remain up during debugging
#endif

	return 0;
}
#endif