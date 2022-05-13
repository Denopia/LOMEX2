#ifndef FILE_READER_H
#define FILE_READER_H

#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <deque>
#include <map>
#include <random>

using namespace std;



class FastqFileReader
{

private:
	std::ifstream fastq_file;
	std::string fastq_line;
	bool reads_left;
	int read_line;
	int first;
	int last;
	int current_read_number;
	//int total_reads;
	std::deque<char> current_read;

public:
    //FastqFileReader(){}
 
    void initialize_me(std::string file_path, int start_position, int end_position);

    int count_reads(std::string file_path);

	void kill_me();

	void roll_to_next_read();

	bool get_reads_left(){return reads_left;}

	std::deque<char> & get_next_read(){return current_read;}

	int get_read_length(){return current_read.size();}

	bool is_empty(){return current_read.empty();}

	char get_next_char(){return current_read.front();}

	char pop_next_char()
	{
		char ret = current_read.front();
		current_read.pop_front();
		return ret;
	}

	void pop_front_char(){current_read.pop_front();}

	int get_current_read_number(){return current_read_number;}

	//int get_total_reads(){return last;}

};

class FastaFileReader
{

private:
	std::ifstream fasta_file;
	int read_number;
	int character_number;
	int line_position;
	std::string current_line;
	bool im_finished;
	int lines_read;

public:
    //FastaFileReader(){}
 
    void initialize_me(std::string file_path);

	void kill_me();

	std::tuple<int, int, char> give_next_read_position_and_character();

	void read_next_line();

	int get_read_number(){return read_number;}

	bool get_im_finished(){return im_finished;}

};


int count_fastq_lines(std::string reads_path);


#endif