
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>

#include <algorithm>
#include <map>
#include <filesystem>
#include <mutex>

#include <tbb/task_group.h>

namespace fs = std::filesystem;


std::map<int, std::vector<double>> load_lbx_single(std::string path)
{
	// measurements per barcode
	std::map<int, std::vector<double>> result;

	std::ifstream ifs(path);
	ifs.imbue(std::locale("C"));

	// skip header
	{
		std::string line;
		std::getline(ifs, line);
	}

	// load the values
	for (;;) {
		int barcode;
		int fi;

#if 1
		ifs >> barcode >> fi;
#else
		std::getline(ifs, line);
		size_t pivot = line.find('\t');

		barcode = std::atoi(line.c_str());
		fi = std::atoi(line.c_str() + pivot);
#endif
		if (ifs) {
			result[barcode].push_back(fi);
		} else {
			break;
		}
	}

	// sort measured values (per barcode)
	for (auto &kvp: result) {
		std::sort(kvp.second.begin(), kvp.second.end());
	}

	return result;
}


std::map<int, std::map<std::string, std::vector<double>>>
//std::map<std::pair<int, std::string>, std::vector<double>>
load_lbx(std::string root, std::string plate)
{
	std::vector<std::string> wells;

	for (const auto &file: fs::directory_iterator(root)) {
		std::string filename = file.path().filename();
		std::string well = filename.substr(plate.size() + 1, 3);

		wells.push_back(well);
	}

	std::map<int, std::map<std::string, std::vector<double>>> result;
	//std::map<std::pair<int, std::string>, std::vector<double>> result;
	std::mutex mtx;

	tbb::task_group tg;

	for (const auto &well: wells) {
		tg.run([&, well]() {
			auto r = load_lbx_single(root + "/" + plate + "_" + well + ".txt");
			std::unique_lock lock(mtx);

			for (auto &kvp: r) {
				result[kvp.first][well] = std::move(kvp.second);
				//result[{kvp.first, well}] = std::move(kvp.second);
			}
		});
	}

	tg.wait();

	return result;
}


void save_gct(std::string path, const std::map<int, std::map<std::string, double>> &dataset)
{
	std::vector<std::string> wells;

	for (const auto &kvp: dataset.begin()->second) {
		wells.push_back(kvp.first);
	}

	std::sort(wells.begin(), wells.end());

	std::ofstream ofs(path);
	ofs.imbue(std::locale("C"));

	ofs << "#1.3\n";
	ofs << dataset.size() << "\t" << wells.size() << "\t0\t0\n";
	ofs << "id";

	for (const auto &well: wells) {
		ofs << "\t" << well;
	}

	ofs << "\n";

	for (auto &kvp: dataset) {
		int barcode = kvp.first;

		ofs << barcode;

		for (const auto &well: wells) {
			ofs << "\t" << kvp.second.at(well);
		}

		ofs << "\n";
	}
}
