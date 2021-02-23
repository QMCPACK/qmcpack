#if COMPILATION_INSTRUCTIONS
c++ -std=c++14 -O3 -Wall -Wextra -Wfatal-errors $0 -o $0x.x -lstdc++fs && $0x.x $@  ; rm -f $0x.x; exit
#endif


#include<iostream>
#include <experimental/filesystem>
#include<string>
#include<regex>

namespace fs = std::experimental::filesystem;
using namespace std::string_literals;
using std::cout;

int main(int argc, char* argv[]){
	std::vector<fs::path> some;
	if(argc == 1 or argv[1] == "all"s){
		const std::regex communicator_related( "communicator_.*\\.cpp" );
		for(fs::path const& p: fs::directory_iterator("./")){
			if(p == "./test.cpp") continue;
		    std::smatch what;
			std::string tmp =p.filename().string();
			if(not regex_match(tmp, what, communicator_related)) continue;
			some.push_back(p);
		}
	}else if(argc==2 and argv[1] == "some"s){
		some = {
			"./communicator_split.cpp", 
			"./communicator_send.cpp",
			"./communicator_reduce.cpp"
		};
	}
	for(fs::path const& p: some){
		cout << p << '\n';
		std::system(("sh " + p.string()).c_str());
	}
}

