/*

	File reading functionality

*/

#include "file_reader.hpp"
#include "fun_kmers.hpp"
#include "file_writer.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <deque>
#include <map>
#include <random>
#include <math.h>
#include <limits>

using namespace std;


/*
	Function to count fastq file lines
*/
int count_fastq_lines(std::string reads_path)
{
	std::ifstream fastq_file;
	fastq_file.open(reads_path, std::ifstream::in);
	std::string fastq_line;
	int number_of_lines = 0;
	while (std::getline(fastq_file, fastq_line, '\n')){number_of_lines+=1;}
	return number_of_lines;
}


/*
	Class for reading fastq files
*/

void FastqFileReader::initialize_me(std::string file_path, int start, int end)
{
	total_reads = count_reads(file_path);
	reads_left = true;
	first = 0;
	last = total_reads;
	read_line = 0;
	current_read_number = 0;
	fastq_line = "";
	fastq_file.open(file_path, std::ifstream::in);
}


/*
	Count the number of reads in a fastq file
*/

int FastqFileReader::count_reads(std::string file_path)
{
	std::ifstream rfile;
	rfile.open(file_path, std::ifstream::in);
	std::string fileline;
	int totalreads = 0;
	int thisline = 0;
	while (std::getline(rfile, fileline, '\n'))
	{
		thisline+=1;
		if (thisline == 2){totalreads+=1;}
		if (thisline == 4){thisline = 0;}
	}
	rfile.close();
	rfile.clear();
	return totalreads;
}




void FastqFileReader::kill_me()
{
	fastq_file.close();
	fastq_file.clear();
}

void FastqFileReader::roll_to_next_read()
{
	// Read fastq file line by line
	bool new_read_gotten = false;
	char current_nucleotide;
	char push_nucleotide;
	current_read.clear();
	while (std::getline(fastq_file, fastq_line, '\n'))
	{
		//std::cout << fastq_line << "\n";
		// Line with the actual read encountered, it is the second line
		if (read_line == 1)
		{
			current_read_number += 1;
			if (current_read_number > last)
			{
				reads_left = false;
				kill_me();
				break;
			}

			if (current_read_number >= first)
			{
				int linelen = fastq_line.length();
				for (int i = 0; i < linelen; i++)
				{
					current_nucleotide = fastq_line.at(i);
					if (current_nucleotide == 'C' || current_nucleotide == 'c'){push_nucleotide = 'C';}
					else if (current_nucleotide == 'A' || current_nucleotide == 'a'){push_nucleotide = 'A';}
					else if (current_nucleotide == 'T' || current_nucleotide == 't'){push_nucleotide = 'T';}
					else if (current_nucleotide == 'G' || current_nucleotide == 'g'){push_nucleotide = 'G';}
					else {push_nucleotide = 'N';}
					current_read.push_back(push_nucleotide);
				}
				new_read_gotten = true;
				reads_left = true;
			}

			if (current_read_number == last)
			{
				reads_left = false;
				kill_me();
			}

		}
		// Update line counters
		read_line += 1;
		// Reached the end of a read
		if (read_line == 4){read_line = 0;}
		// Break after a new read has been found
		if (new_read_gotten){break;}
	}
	if (!new_read_gotten)
	{
		reads_left = false;
		kill_me();
	}
}






/*
	Class for reading fasta files
*/

void FastaFileReader::initialize_me(std::string file_path)
{
	read_number = 0;
	character_number = 0;
	line_position = 0;
	current_line = "\n";
	im_finished = false;
	lines_read = 0;
	fasta_file.open(file_path, std::ifstream::in);
}

void FastaFileReader::kill_me()
{
	fasta_file.close();
	fasta_file.clear();
}


std::tuple<int, int, char> FastaFileReader::give_next_read_position_and_character()
{

	if (line_position > current_line.length()-1)
	{
		read_next_line();
	}
	if (line_position < current_line.length())
	{
		while((current_line.at(line_position) == '\n') || (current_line.at(line_position) == '>'))
		{
			if (current_line.at(line_position) == '>')
			{
				read_number += 1;
				character_number = 0;
			}
			if (im_finished){break;}
			read_next_line();
		}
	}

	char return_char;
	int return_read;
	int return_position;

	if (im_finished)
	{
		kill_me();
		return_char = 'Z';
		return_read = -1;
		return_position = -1;
	} 
	else
	{
		return_char = current_line.at(line_position);
		return_read = read_number;
		return_position = character_number;
	}

	line_position += 1;
	character_number += 1;

	return std::make_tuple(return_read, return_position, return_char);

}


void FastaFileReader::read_next_line()
{
	std::getline(fasta_file, current_line, '\n');
	//std::cout << current_line << std::endl;
	line_position = 0;
	lines_read += 1;
	if (current_line == "")
	{
		im_finished = true;
	}
}




