#include "fun_poa.hpp"
#include <iostream>
#include <vector>
#include <stdio.h>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <chrono>

/*extern "C" 
{
 	#include "lpo.h"
	#include "msa_format.h"
	#include "align_score.h"
}*/


int xlate_lpo_to_al(LPOSequence_T *seq,
		    int nsymbol,char symbol[],int ibundle,
		    char gap_character,
		    char ***p_seq_pos,char **p_p,char **p_include)
{
  int i,j,iring=0,nring=0,current_ring=0,iprint;
  char **seq_pos=NULL,*p=NULL,*include_in_save=NULL;
  LPOLetterSource_T *source;

  LOOPF (i,seq->length) /* COUNT TOTAL #ALIGNMENT RINGS IN THE LPO */
    if (seq->letter[i].ring_id != current_ring) { /* NEXT RING */
      current_ring=seq->letter[i].ring_id;
      nring++;
    }
  nring++; /* DON'T FORGET TO COUNT THE LAST RING!!! */
  
  CALLOC(seq_pos,seq->nsource_seq,char *); /* ALLOCATE MAP ARRAY*/
  CALLOC(p,seq->nsource_seq*nring,char);
  LOOP (i,seq->nsource_seq) /* BUILD POINTER ARRAY INTO MAP ARRAY */
    seq_pos[i]=p+i*nring;
  memset(p,gap_character,seq->nsource_seq*nring);
  /* DEFAULT IS NO SEQUENCE PRESENT AT THIS POSITION */

  current_ring=0; /* RESET TO BEGINNING */
  LOOPF (i,seq->length) { /* NOW MAP THE LPO TO A FLAT LINEAR ORDER */
    if (seq->letter[i].ring_id != current_ring) { /* NEXT RING */
      current_ring=seq->letter[i].ring_id;
      iring++;
    } /* MAP EACH SOURCE SEQ ONTO LINEAR ORDER INDEXED BY iring */
    for (source= &seq->letter[i].source;source;source=source->more)
      if (symbol && seq->letter[i].letter<nsymbol)  /* TRANSLATE TO symbol */
	seq_pos[source->iseq][iring]= symbol[seq->letter[i].letter];
      else  /* NO NEED TO TRANSLATE */
	seq_pos[source->iseq][iring]= seq->letter[i].letter;
  }

  if (ibundle>=0) { /* ONLY SAVE SEQS THAT ARE IN THIS BUNDLE */
    CALLOC(include_in_save,nring,char); /* BLANK FLAGS: WHAT RINGS TO SHOW*/
    LOOP (iring,nring) { /* CHECK EACH RING TO SEE IF IT'S IN BUNDLE */
      LOOP (i,seq->nsource_seq) {
	if (seq_pos[i][iring]!=gap_character /* ALIGNED HERE! */
	    && seq->source_seq[i].bundle_id == ibundle) { /* PART OF BUNDLE!*/
	  include_in_save[iring]=1; /* SO INCLUDE THIS RING */
	  break;
	}
      }
    }
  }
  //if (p_seq_pos)
  //  *p_seq_pos = seq_pos;
  //else
  //  FREE (seq_pos);
  //FREE (p);
  //FREE (include_in_save);

   if (p_seq_pos)
    *p_seq_pos = seq_pos;
  else
    FREE (seq_pos);
  
  if (p_p)
    *p_p = p;
  else
    FREE (p);
  
  if (p_include)
    *p_include = include_in_save;
  else
    FREE (include_in_save);

  //FREE()
  return nring;
}

//#define SEQ_NAME_MAX 50
//#define SEQ_TITLE_MAX 50
//#define SEQ_LENGTH_MAX 100

void LOMEXPOAGlue::initialize_me(std::string score_matrix_path)
{
	//std::cout << "INIT START\n";
	int path_len = score_matrix_path.length();
	char score_matrix_path_ca[path_len+1];
	strcpy(score_matrix_path_ca, score_matrix_path.c_str());
	score_matrix_int = read_score_matrix(score_matrix_path_ca, &score_matrix);
	//std::cout << "INIT END\n";
}


int LOMEXPOAGlue::poa_consensus_sequences(Sequence_T **mainseqs, std::vector<std::string>& input_sequences, float bundling_threshold, float support_threshold, int min_bundle_size_threshold, int min_char_support_threshold)
{

	int seq_name_max = 50;
	int seq_title_max = 50;
	int seq_length_max = 100;

	char seq_name[seq_name_max]="";
	char seq_line[seq_length_max]="";
	char seq_title[seq_title_max]="";

	stringptr tmp_seq=STRINGPTR_EMPTY_INIT;

	int length = 0;
	int nseq = 0;

	//LPOSequence_T *seqC = NULL;

	//std::cout << "GLUE 1\n";
	consensus_sequences.clear();

	//LPOSequence_T **input_seqs = NULL;
	//LPOSequence_T *seq = NULL;
	//LPOSequence_T *lpo_out = NULL;
	int n_input_seqs = 0;
	int max_input_seqs = 50;
	int isid = 0;
	int inseqnum = input_sequences.size();
	
	// Options ??
	int do_switch_case = 1; /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WHAT SHOULD I PUT HERE??
	int use_aggressive_fusion = 0;
	int do_progressive = 0;
	int do_global = 1;
	int do_preserve_sequence_order = 0;
	//char * pair_score_file = NULL;

	std::string name_str;
	std::string title_str;
	
	//char name[100];

	//std::cout << "GLUE 1\n";

	// Read input sequences
	for (int iis = 0; iis < input_sequences.size(); iis+=1)
	//for (auto seq_str : input_sequences)
	{
    	name_str = "tn" + std::to_string(isid);
    	title_str = "tt" + std::to_string(isid);
    	//std::cout << "GLUE 2\n";
    	//seq_name = &name_str[0];
		//seq_title = &title_str[0];
		//seq_line = &(input_sequences[iis])[0];
		strcpy(seq_name, name_str.c_str());
		strcpy(seq_title, title_str.c_str());
		strcpy(seq_line, input_sequences[iis].c_str());
		//std::cout << "GLUE 3\n";
		stringptr_cat_pos(&tmp_seq,seq_line,&length);
		//std::cout << "GLUE 4\n";
		//std::cout << "ORIGINAL " << seq_line << "\n";
		//std::cout << "TRANSFORMED " << tmp_seq.p << "\n";
		//std::cout << "NAME:" << test_name << " TITLE: " << test_title << " SEQ: " << new_seq << "\n"; 
		//create_seq(isid, &seq, name, title, sequ, do_switch_case); //new
		if (create_seq(nseq,mainseqs,seq_name,seq_title,tmp_seq.p,do_switch_case)){nseq++;}
		isid += 1;
		//std::cout << "GLUE 5\n";
		//nseq += 1;
	}



	//std::cout << "GLUE 6\n";
	stringptr_free(&tmp_seq);
	
	//CALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
	//FREE(seq);
	//FREE (input_seqs);
	//if (nseq > 0){FREE (seq);}

	//if (nseq > 0){FREE (seqC);}

	return nseq;

	/*

	//std::cout << "GLUE 3\n";


	//std::cout << "START PUTTING SEQUENCES INTO LPO\n";

	// Make LPO? out of the input sequences
	for (int i = 0; i < nseq; i+=1) {
		input_seqs[n_input_seqs++] = &(seq[i]);
		initialize_seqs_as_lpo(1,&(seq[i]),&score_matrix);
		if (n_input_seqs == max_input_seqs) {
			max_input_seqs *= 2;
			REALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
		}
	}


	//FREE (input_seqs);
	//if (nseq > 0){FREE (seq);}

	//return;





	// Bundle consensus sequences
	lpo_out = buildup_progressive_lpo (isid, input_seqs, &score_matrix, 
		use_aggressive_fusion, do_progressive, pair_score_file, 
		matrix_scoring_function, do_global, do_preserve_sequence_order);

	FREE(lpo_out->title);
	
	generate_lpo_bundles(lpo_out, bundling_threshold);

	char gap_char = '-';
	int nring = 0;
	int ibundle = 1;
	// What are these?
	char **seq_pos = NULL;
	char *p = NULL;
	char *include_in_save = NULL;
	
	// xlate means ???
	//std::cout << "START DOING NRING CALCULATIONS\n";


	nring = xlate_lpo_to_al(lpo_out, score_matrix.nsymbol, score_matrix.symbol, ibundle, gap_char, &seq_pos, &p, &include_in_save);


 	//nring = xlate_lpo_to_al(seq, nsymbol, symbol, ALL_BUNDLES, gap_char, &seq_pos, &p, &include_in_save);
	//FREE(p); / DUMP TEMPORARY MEMORY /
	//FREE(include_in_save);
	//FREE(seq_pos);




	//std::cout << "GLUE 7\n";
	//std::cout << "NRING THING DONE\n";
	// Extract sequences
	//std::cout << "****************** CONSENSUS SEQUENCES ************************\n";
	//for (int i = inseqnum; i < lpo_out->nsource_seq; i += 1)

	std::string conseq;
	std::unordered_map<int,int> bundle_counts;
	std::vector<std::string> input_seqs_aligned;
	std::unordered_map<int, std::vector<std::string>> input_seq_bundles;
	std::vector<std::string> consensus_seqs_aligned;


	for (int i = 0; i < lpo_out->nsource_seq; i += 1)
	{
		conseq = "";
		//std::cout << "*****\n";
		//std::cout << "SEQ NAME: " << lpo_out->source_seq[i].name << "\n";
		//std::cout << "SEQ TITLE: " << lpo_out->source_seq[i].title << "\n";
		//std::cout << "SEQ BUNDLE: " << lpo_out->source_seq[i].bundle_id << "\n";
		
		// Get sequence's bundle
		int seq_bundle = lpo_out->source_seq[i].bundle_id;

		// Build a string for the aligned sequence
		for (int j = 0; j < nring; j+=1){conseq = conseq + seq_pos[i][j];}
		
		// Input sequence bundles
		if (i < inseqnum)
		{
			if (seq_bundle >= 0)
			{
				if (bundle_counts.count(seq_bundle) == 0){
					bundle_counts[seq_bundle] = 0;
				}
				bundle_counts[seq_bundle] += 1;
				input_seq_bundles[seq_bundle].push_back(conseq);
			}		
			//continue;
		}

		// Print debug stuff
		/*
		if (i < inseqnum){
			if (lpo_out->source_seq[i].bundle_id >= 0){
				std::cout << lpo_out->source_seq[i].bundle_id << "  " << conseq << "\n";
			} else {
				std::cout << lpo_out->source_seq[i].bundle_id << " " << conseq << "\n";	
			}
		} 
		*/

		/*
		
		// Determine if consensus sequence is good enough to be put into the 
		
		// This must be a consensus sequence
		if (i >= inseqnum)
		{
			bool is_real_bundle = bundle_counts.count(seq_bundle) > 0;
			bool is_supported_percent = bundle_counts[seq_bundle] > support_threshold*inseqnum;
			bool is_supported_hard = bundle_counts[seq_bundle] > min_bundle_size_threshold;

			bool enough_support = true;
			std::vector<int> character_supports;

			for (int n = 0; n < nring; n+=1)
			{
				int position_support = 0;

				// Deletions (dashses) do not need to be supported
				if (conseq.at(n) == '-')
				{
					character_supports.push_back(0);
					continue;
				}

				for (int m = 0; m < input_seq_bundles[seq_bundle].size(); m+=1)
				{
					if (input_seq_bundles[seq_bundle][m].at(n) == conseq.at(n)){position_support += 1;}
				}

				if (position_support < min_char_support_threshold)
				{
					character_supports.push_back(-1);
					//enough_support = false;
					//break;
				} else {
					character_supports.push_back(1);
				}
			}


			// Find supported start
			int supported_start;
			for (int csi = 0; csi < nring; csi+=1)
			{
				if (character_supports.at(csi) == 1){supported_start = csi;break;}
			}

			// Find supported end
			int supported_end;
			for (int csi = nring-1; csi >= 0; csi-=1)
			{
				if (character_supports.at(csi) == 1){supported_end = csi;break;}
			}
			
			// No random unsupported characters in the middle
			int from_support_to_no_support = 0;
			int last_non_deletion = character_supports.at(0);
			for (int csi = 1; csi < nring; csi+=1)
			{
				if ((character_supports.at(csi) == -1) && (last_non_deletion == 1)){from_support_to_no_support+=1;}
				if (character_supports.at(csi) == -1){last_non_deletion = -1;}
				if (character_supports.at(csi) == 1){last_non_deletion = 1;}
			}

			// Determine if enough support
			if (from_support_to_no_support > 1){enough_support = false;}

			// Print accordingly
			if (enough_support && is_real_bundle && is_supported_percent && is_supported_hard){
				consensus_sequences.push_back(conseq);	
				//std::cout << lpo_out->source_seq[i].bundle_id << "  " << conseq << " ** ("<< supported_start << "," << supported_end << ")\n";
			} else {
				continue;
				//std::cout << lpo_out->source_seq[i].bundle_id << "  " << conseq << " * \n";
			}

			character_supports.clear();

		}

	}

	//std::cout << "GLUE 8\n";

	// Free allocated memory

	//std::cout << "FREEING ALLOCATED MEMORY, ERROR MIGHTS OCCUR!!!\n";

	int ii;
	int jj;
	for (ii = 0; ii < n_input_seqs; ii+=1) {
		for (jj = 0; jj < nseq; jj+=1) {
			if (input_seqs[ii]==&(seq[jj])) {
				break;
			}
		}
		free_lpo_sequence(input_seqs[ii],(jj==nseq));
	}


	//std::cout << "FREEING MORE ALLOCATED MEMORY, ERROR MIGHTS OCCUR!!!\n";


	FREE (input_seqs);
	if (nseq > 0){FREE (seq);}
	//FREE (lpo_out); /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DO I NEED THIS?

	consensus_seqs_aligned.clear();
	input_seq_bundles.clear();
	input_seqs_aligned.clear();
	bundle_counts.clear();


	FREE(p);
	FREE(include_in_save);
	FREE(seq_pos);
	
	// Return consensus sequences
	//return consensus_sequences;

	*/


}

// THE MAIN FUNCTION THAT IS USED!!!!!!!!!!!!!!!!!!!!!!!!!!
std::tuple<int, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> construct_consensus_seqs(std::vector<std::string>& occurrence_sequences, std::vector<std::string>& consensus_sequences, ResidueScoreMatrix_T *score_matrix, int do_switch_case,
	float bundling_threshold, float support_threshold, int min_bundle_size_threshold, int min_char_support_threshold)
{

	std::chrono::high_resolution_clock::time_point poa_start_time, poa_end_time;
	std::chrono::high_resolution_clock::time_point poa_build_start_time, poa_build_end_time;
	std::chrono::high_resolution_clock::time_point poa_cons_start_time, poa_cons_end_time;

	std::chrono::high_resolution_clock::time_point poa_p1_start_time, poa_p1_end_time;
	std::chrono::high_resolution_clock::time_point poa_p2_start_time, poa_p2_end_time;
	std::chrono::high_resolution_clock::time_point poa_p3_start_time, poa_p3_end_time;


	std::chrono::high_resolution_clock::time_point filtering_start_time, filtering_end_time;

	poa_start_time = std::chrono::high_resolution_clock::now();
	poa_build_start_time = std::chrono::high_resolution_clock::now();

	LPOSequence_T *mainseq = NULL;

	int seq_name_max = 50;
	int seq_title_max = 50;
	int seq_length_max = 100;

	char seq_name[seq_name_max]="";
	char seq_line[seq_length_max]="";
	char seq_title[seq_title_max]="";

	stringptr tmp_seq=STRINGPTR_EMPTY_INIT;

	int length = 0;
	int nseq = 0;
	int isid = 0;
	int inseqnum = occurrence_sequences.size();
	
	std::string name_str = "";
	std::string title_str = "";
	std::string line_str = "";
	

	poa_p1_start_time = std::chrono::high_resolution_clock::now();

	// Read input sequences
	for (int iis = 0; iis < occurrence_sequences.size(); iis+=1)
	{
    	name_str = "tn" + std::to_string(isid);
    	title_str = "tt" + std::to_string(isid);
		//line_str = occurrence_sequences[iis];
		
		strcpy(seq_name, name_str.c_str());
		strcpy(seq_title, title_str.c_str());
		strcpy(seq_line, occurrence_sequences[iis].c_str());

		stringptr_cat_pos(&tmp_seq, seq_line, &length);
		if (create_seq(nseq, &mainseq, seq_name, seq_title, tmp_seq.p, do_switch_case)){nseq++;} ///   THIS IS THE PROBLEM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! tmp_sq.p???
		//if (create_seq(nseq, &mainseq, seq_name, seq_title, seq_line, do_switch_case)){nseq+=1;}
		
		if (tmp_seq.p){
			tmp_seq.p[0]='\0';
			length = 0;
		}
		isid += 1;
	}

	stringptr_free(&tmp_seq);

	poa_p1_end_time = std::chrono::high_resolution_clock::now();


	poa_p2_start_time = std::chrono::high_resolution_clock::now();


	// Use strings to create lpo
	LPOSequence_T **input_seqs = NULL;
	
	int n_input_seqs = 0;
	int max_input_seqs = 40;  //                                                       WAS 10!!!!!!!!!!!!!!

	CALLOC (input_seqs, max_input_seqs, LPOSequence_T *);

	int di = 0;
	for (di=0; di<nseq; di++) {
		input_seqs[n_input_seqs++] = &(mainseq[di]);
		initialize_seqs_as_lpo(1,&(mainseq[di]),score_matrix);
		if (n_input_seqs == max_input_seqs) {
			max_input_seqs *= 2;
			REALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
		}
	}

	poa_p2_end_time = std::chrono::high_resolution_clock::now();


	// Bottleneck
	//====================================================================================
	poa_p3_start_time = std::chrono::high_resolution_clock::now();

	int use_aggressive_fusion = 0;
	int do_progressive = 0;
	int do_global = 1;
	int do_preserve_sequence_order = 0;
	char * pair_score_file = NULL;
	//float bundling_threshold = 0.4;

	LPOSequence_T *lpo_out = NULL;

	// Build progressive lpo
	lpo_out = buildup_progressive_lpo (isid, input_seqs, score_matrix, 
		use_aggressive_fusion, do_progressive, pair_score_file, 
		matrix_scoring_function, do_global, do_preserve_sequence_order);

	poa_p3_end_time = std::chrono::high_resolution_clock::now();
	//====================================================================================



	poa_build_end_time = std::chrono::high_resolution_clock::now();

	poa_cons_start_time = std::chrono::high_resolution_clock::now();

	//FREE (lpo_out->title);

	// Create lpo bundles
	generate_lpo_bundles(lpo_out, bundling_threshold);

	char gap_char = '-';
	int nring;
	int ibundle = -1;
	char **seq_pos = NULL;
	char *p = NULL;
	char *include_in_save = NULL;
	

	nring = xlate_lpo_to_al(lpo_out, score_matrix->nsymbol, score_matrix->symbol, ALL_BUNDLES, gap_char, &seq_pos, &p, &include_in_save);

	


/*   THIS BLOCK WORKS NICELY

	std::string conseq;

	// This is the simple way to pick everything
	for (int i = 0; i < lpo_out->nsource_seq; i += 1)
	{
		conseq = "";
		int seq_bundle = lpo_out->source_seq[i].bundle_id;
		if (i < inseqnum){continue;}
		for (int j = 0; j < nring; j+=1){conseq = conseq + seq_pos[i][j];}
		consensus_sequences.push_back(conseq);
	}

	//std::cout << "DELETING INPUTS 1\n";


	// FREE ALLOCATED MEMORY
	int mi=0;
	int mj=0;
	for (mi=0;mi<n_input_seqs;mi++) {
		for (mj=0;mj<nseq;mj++) {
			if (input_seqs[mi]==&(mainseq[mj])){break;}
		}
		//free_lpo_sequence(input_seqs[mi],(mj==nseq));
		//free_lpo_sequence(mainseqs[mi],(mj==nseq));
		//std::cout << "FREE 1\n";
		free_lpo_sequence(input_seqs[mi],(mj==nseq));
		//std::cout << "FREE 2\n";
		//free_lpo_sequence(&(mainseq[mi]),false);
	}

	//std::cout << "DELETING INPUTS 2\n";
	FREE (input_seqs);
	//std::cout << "DELETING OUTS\n";
	
	//free_lpo_sequence(lpo_out, true);                        IS IT THIS OR BELOW?
	if (nseq > 0){ FREE (mainseq);}
	//FREE (lpo_out);
	//std::cout << "NOT DELETING MAIN\n";
	//if (nseq>0) FREE (mainseq);

	// DUMP TEMPORARY MEMORY
	FREE(p); 
  FREE(include_in_save);
  FREE(seq_pos);

	return nseq;
*/

	poa_cons_end_time = std::chrono::high_resolution_clock::now();

	filtering_start_time = std::chrono::high_resolution_clock::now();

	std::unordered_map<int,int> bundle_counts;
	std::vector<std::string> input_seqs_aligned;
	std::unordered_map<int, std::vector<std::string>> input_seq_bundles;
	std::vector<std::string> consensus_seqs_aligned;
	std::string conseq;


	for (int i = 0; i < lpo_out->nsource_seq; i += 1)
	{
		conseq = "";
		//std::cout << "*****\n";
		//std::cout << "SEQ NAME: " << lpo_out->source_seq[i].name << "\n";
		//std::cout << "SEQ TITLE: " << lpo_out->source_seq[i].title << "\n";
		//std::cout << "SEQ BUNDLE: " << lpo_out->source_seq[i].bundle_id << "\n";
		
		// Get sequence's bundle
		int seq_bundle = lpo_out->source_seq[i].bundle_id;

		// Build a string for the aligned sequence
		for (int j = 0; j < nring; j+=1){conseq = conseq + seq_pos[i][j];}
		
		// Input sequence bundles
		if (i < inseqnum)
		{
			if (seq_bundle >= 0)
			{
				if (bundle_counts.count(seq_bundle) == 0){
					bundle_counts[seq_bundle] = 0;
				}
				bundle_counts[seq_bundle] += 1;
				input_seq_bundles[seq_bundle].push_back(conseq);
			}		
			//continue;
		}

		// Print debug stuff
		/*
		if (i < inseqnum){
			if (lpo_out->source_seq[i].bundle_id >= 0){
				std::cout << lpo_out->source_seq[i].bundle_id << "  " << conseq << "\n";
			} else {
				std::cout << lpo_out->source_seq[i].bundle_id << " " << conseq << "\n";	
			}
		} 
		*/

		// Determine if consensus sequence is good enough to be put into the 
		// This must be a consensus sequence

		if (i >= inseqnum)
		{
			bool is_real_bundle = bundle_counts.count(seq_bundle) > 0;
			bool is_supported_percent = bundle_counts[seq_bundle] > support_threshold*inseqnum;
			bool is_supported_hard = bundle_counts[seq_bundle] > min_bundle_size_threshold;

			bool enough_support = true;
			std::vector<int> character_supports;

			for (int n = 0; n < nring; n+=1)
			{
				int position_support = 0;

				// Deletions (dashses) do not need to be supported
				if (conseq.at(n) == '-')
				{
					character_supports.push_back(0);
					continue;
				}

				for (int m = 0; m < input_seq_bundles[seq_bundle].size(); m+=1)
				{
					if (input_seq_bundles[seq_bundle][m].at(n) == conseq.at(n)){position_support += 1;}
				}

				if (position_support < min_char_support_threshold)
				{
					character_supports.push_back(-1);
					//enough_support = false;
					//break;
				} else {
					character_supports.push_back(1);
				}
			}


			// Find supported start
			int supported_start;
			for (int csi = 0; csi < nring; csi+=1)
			{
				if (character_supports.at(csi) == 1){supported_start = csi;break;}
			}

			// Find supported end
			int supported_end;
			for (int csi = nring-1; csi >= 0; csi-=1)
			{
				if (character_supports.at(csi) == 1){supported_end = csi;break;}
			}
			
			// No random unsupported characters in the middle
			int from_support_to_no_support = 0;
			int last_non_deletion = character_supports.at(0);
			for (int csi = 1; csi < nring; csi+=1)
			{
				if ((character_supports.at(csi) == -1) && (last_non_deletion == 1)){from_support_to_no_support+=1;}
				if (character_supports.at(csi) == -1){last_non_deletion = -1;}
				if (character_supports.at(csi) == 1){last_non_deletion = 1;}
			}

			// Determine if enough support
			if (from_support_to_no_support > 1){enough_support = false;}

			// Print accordingly
			if (enough_support && is_real_bundle && is_supported_percent && is_supported_hard){
				consensus_sequences.push_back(conseq);	
				//std::cout << lpo_out->source_seq[i].bundle_id << "  " << conseq << " ** ("<< supported_start << "," << supported_end << ")\n";
			} else {
				continue;
				//std::cout << lpo_out->source_seq[i].bundle_id << "  " << conseq << " * \n";
			}

			character_supports.clear();
		}
	}

	filtering_end_time = std::chrono::high_resolution_clock::now();


	// FREE ALLOCATED MEMORY
	int mi=0;
	int mj=0;
	for (mi=0;mi<n_input_seqs;mi++) {
		for (mj=0;mj<nseq;mj++) {
			if (input_seqs[mi]==&(mainseq[mj])){break;}
		}
		//free_lpo_sequence(input_seqs[mi],(mj==nseq));
		//free_lpo_sequence(mainseqs[mi],(mj==nseq));
		//std::cout << "FREE 1\n";
		free_lpo_sequence(input_seqs[mi],(mj==nseq));
		//std::cout << "FREE 2\n";
		//free_lpo_sequence(&(mainseq[mi]),false);
	}

	//std::cout << "DELETING INPUTS 2\n";
	FREE (input_seqs);
	//std::cout << "DELETING OUTS\n";
	
	//free_lpo_sequence(lpo_out, true);                        IS IT THIS OR BELOW?
	if (nseq > 0){ FREE (mainseq);}
	//FREE (lpo_out);
	//std::cout << "NOT DELETING MAIN\n";
	//if (nseq>0) FREE (mainseq);

	// DUMP TEMPORARY MEMORY
	FREE(p); 
  FREE(include_in_save);
  FREE(seq_pos);


  poa_end_time = std::chrono::high_resolution_clock::now();

  uint64_t poa_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(poa_end_time-poa_start_time).count();
  uint64_t filtering_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(filtering_end_time-filtering_start_time).count();

  uint64_t poa_build_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(poa_build_end_time-poa_build_start_time).count();
  uint64_t poa_cons_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(poa_cons_end_time-poa_cons_start_time).count();

  int true_poa_microseconds = poa_microseconds - filtering_microseconds;

  uint64_t poa_p1_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(poa_p1_end_time-poa_p1_start_time).count();
  uint64_t poa_p2_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(poa_p2_end_time-poa_p2_start_time).count();
  uint64_t poa_p3_microseconds = std::chrono::duration_cast<std::chrono::microseconds>(poa_p3_end_time-poa_p3_start_time).count();


	return std::make_tuple(nseq, true_poa_microseconds, filtering_microseconds, poa_build_microseconds, poa_cons_microseconds, poa_p1_microseconds, poa_p2_microseconds, poa_p3_microseconds);

	/* Works from here */
	// FREE ALLOCATED MEMORY
	//int mi=0;
	//int mj=0;
	//for (mi=0;mi<n_input_seqs;mi++) {
	//	for (mj=0;mj<nseq;mj++) {
	//		if (input_seqs[mi]==&(mainseq[mj])){break;}
	//	}
	//	free_lpo_sequence(input_seqs[mi],(mj==nseq));
	//}
	//FREE (input_seqs);
	//if (nseq>0) FREE (mainseq);
	//return nseq;
}




void LOMEXPOAGlue::poa_consensus_sequences_unrefined(std::vector<std::string>& input_sequences, float bundling_threshold, float support_threshold, int min_bundle_size_threshold, int min_char_support_threshold)
{

	//std::cout << "GLUE 1\n";
	consensus_sequences.clear();

	LPOSequence_T **input_seqs = NULL;
	LPOSequence_T *seq = NULL;
	LPOSequence_T *lpo_out = NULL;
	int n_input_seqs = 0;
	int max_input_seqs = 50;
	int nseq = 0;
	int isid = 0;
	int inseqnum = input_sequences.size();
	
	// Options ??
	int do_switch_case = 1; /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WHAT SHOULD I PUT HERE??
	int use_aggressive_fusion = 0;
	int do_progressive = 0;
	int do_global = 1;
	int do_preserve_sequence_order = 0;
	char * pair_score_file = NULL;

	CALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
			
	char * name;
	char * title;
	char * sequ;

	//std::cout << "GLUE 2\n";

	// Read input sequences
	for (auto seq_str : input_sequences)
	{
    	std::string name_str("tn"+std::to_string(isid));
    	std::string title_str("tt"+std::to_string(isid));
		name = &name_str[0];
		title = &title_str[0];
		sequ = &seq_str[0];
		//std::cout << "NAME:" << test_name << " TITLE: " << test_title << " SEQ: " << new_seq << "\n"; 
		create_seq(isid, &seq, name, title, sequ, do_switch_case);
		isid += 1;
		nseq += 1;
	}

	return;

	//std::cout << "GLUE 3\n";


	//std::cout << "START PUTTING SEQUENCES INTO LPO\n";

	// Make LPO? out of the input sequences
	for (int i = 0; i < nseq; i+=1) {
		input_seqs[n_input_seqs++] = &(seq[i]);
		initialize_seqs_as_lpo(1,&(seq[i]),&score_matrix);
		if (n_input_seqs == max_input_seqs) {
			max_input_seqs *= 2;
			REALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
		}
	}

	// Bundle consensus sequences
	lpo_out = buildup_progressive_lpo (isid, input_seqs, &score_matrix, 
		use_aggressive_fusion, do_progressive, pair_score_file, 
		matrix_scoring_function, do_global, do_preserve_sequence_order);

	FREE (lpo_out->title);
	
	generate_lpo_bundles(lpo_out, bundling_threshold);

	char gap_char = '-';
	int nring;
	int ibundle = -1;
	// What are these?
	char **seq_pos = NULL;
	char *p = NULL;
	char *include_in_save = NULL;
	
	// xlate means ???
	//std::cout << "START DOING NRING CALCULATIONS\n";


	//nring = xlate_lpo_to_al(lpo_out, score_matrix.nsymbol, score_matrix.symbol, ibundle, gap_char, &seq_pos, &p, &include_in_save);

	nring = xlate_lpo_to_al(lpo_out, score_matrix.nsymbol, score_matrix.symbol, ALL_BUNDLES, gap_char, &seq_pos, &p, &include_in_save);


	//std::cout << "GLUE 7\n";
	//std::cout << "NRING THING DONE\n";
	// Extract sequences
	//std::cout << "****************** CONSENSUS SEQUENCES ************************\n";
	//for (int i = inseqnum; i < lpo_out->nsource_seq; i += 1)

	std::string conseq;
	
	/*
	std::unordered_map<int,int> bundle_counts;
	std::vector<std::string> input_seqs_aligned;
	std::unordered_map<int, std::vector<std::string>> input_seq_bundles;
	std::vector<std::string> consensus_seqs_aligned;
	*/


	for (int i = 0; i < lpo_out->nsource_seq; i += 1)
	{
		conseq = "";
		//std::cout << "*****\n";
		//std::cout << "SEQ NAME: " << lpo_out->source_seq[i].name << "\n";
		//std::cout << "SEQ TITLE: " << lpo_out->source_seq[i].title << "\n";
		//std::cout << "SEQ BUNDLE: " << lpo_out->source_seq[i].bundle_id << "\n";
		
		// Get sequence's bundle
		int seq_bundle = lpo_out->source_seq[i].bundle_id;

		// Build a string for the aligned sequence
		
		// Input sequence bundles
		if (i < inseqnum)
		{
			/*
			if (seq_bundle >= 0)
			{
				if (bundle_counts.count(seq_bundle) == 0){
					bundle_counts[seq_bundle] = 0;
				}
				bundle_counts[seq_bundle] += 1;
				input_seq_bundles[seq_bundle].push_back(conseq);
			}
			*/		
			continue;
		}

		for (int j = 0; j < nring; j+=1){conseq = conseq + seq_pos[i][j];}

		consensus_sequences.push_back(conseq);

		// Print debug stuff
		/*
		if (i < inseqnum){
			if (lpo_out->source_seq[i].bundle_id >= 0){
				std::cout << lpo_out->source_seq[i].bundle_id << "  " << conseq << "\n";
			} else {
				std::cout << lpo_out->source_seq[i].bundle_id << " " << conseq << "\n";	
			}
		} 
		*/
		
		// Determine if consensus sequence is good enough to be put into the 
		
		// This must be a consensus sequence
		
		/*
		if (i >= inseqnum)
		{
			bool is_real_bundle = bundle_counts.count(seq_bundle) > 0;
			bool is_supported_percent = bundle_counts[seq_bundle] > support_threshold*inseqnum;
			bool is_supported_hard = bundle_counts[seq_bundle] > min_bundle_size_threshold;

			bool enough_support = true;
			std::vector<int> character_supports;

			for (int n = 0; n < nring; n+=1)
			{
				int position_support = 0;

				// Deletions (dashses) do not need to be supported
				if (conseq.at(n) == '-')
				{
					character_supports.push_back(0);
					continue;
				}

				for (int m = 0; m < input_seq_bundles[seq_bundle].size(); m+=1)
				{
					if (input_seq_bundles[seq_bundle][m].at(n) == conseq.at(n)){position_support += 1;}
				}

				if (position_support < min_char_support_threshold)
				{
					character_supports.push_back(-1);
					//enough_support = false;
					//break;
				} else {
					character_supports.push_back(1);
				}
			}


			// Find supported start
			int supported_start;
			for (int csi = 0; csi < nring; csi+=1)
			{
				if (character_supports.at(csi) == 1){supported_start = csi;break;}
			}

			// Find supported end
			int supported_end;
			for (int csi = nring-1; csi >= 0; csi-=1)
			{
				if (character_supports.at(csi) == 1){supported_end = csi;break;}
			}
			
			// No random unsupported characters in the middle
			int from_support_to_no_support = 0;
			int last_non_deletion = character_supports.at(0);
			for (int csi = 1; csi < nring; csi+=1)
			{
				if ((character_supports.at(csi) == -1) && (last_non_deletion == 1)){from_support_to_no_support+=1;}
				if (character_supports.at(csi) == -1){last_non_deletion = -1;}
				if (character_supports.at(csi) == 1){last_non_deletion = 1;}
			}

			// Determine if enough support
			if (from_support_to_no_support > 1){enough_support = false;}

			// Print accordingly
			if (enough_support && is_real_bundle && is_supported_percent && is_supported_hard){
				consensus_sequences.push_back(conseq);	
				//std::cout << lpo_out->source_seq[i].bundle_id << "  " << conseq << " ** ("<< supported_start << "," << supported_end << ")\n";
			} else {
				continue;
				//std::cout << lpo_out->source_seq[i].bundle_id << "  " << conseq << " * \n";
			}

			character_supports.clear();

		}
		*/

	}

	//std::cout << "GLUE 8\n";

	// Free allocated memory

	//std::cout << "FREEING ALLOCATED MEMORY, ERROR MIGHTS OCCUR!!!\n";

	int ii;
	int jj;
	for (ii = 0; ii < n_input_seqs; ii+=1) {
		for (jj = 0; jj < nseq; jj+=1) {
			if (input_seqs[ii]==&(seq[jj])) {
				break;
			}
		}
		free_lpo_sequence(input_seqs[ii],(jj==nseq));
	}

	//std::cout << "FREEING MORE ALLOCATED MEMORY, ERROR MIGHTS OCCUR!!!\n";

	FREE (input_seqs);
	if (nseq > 0){FREE (seq);}



	FREE (p);
	FREE (include_in_save);
	FREE (seq_pos);


	//FREE (lpo_out); /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DO I NEED THIS?

	//consensus_seqs_aligned.clear();
	//input_seq_bundles.clear();
	//input_seqs_aligned.clear();
	//bundle_counts.clear();
	
	// Return consensus sequences
	//return consensus_sequences;

	return;
}


void do_poa_stuff()
{

	char fake_path[] = "../poaV2FREESH/blosum80.mat";

	ResidueScoreMatrix_T score_matrix;
	
	int a = read_score_matrix(fake_path, &score_matrix);

	char buf[1024]; 
	sprintf(buf, "X-Gap Penalties (Open, Aff1, Aff2; LTrunc, LDecay): %d %d %d %d %d \n",
	    score_matrix.gap_penalty_set[0][0],
	    score_matrix.gap_penalty_set[0][1],
	    score_matrix.gap_penalty_set[0][2],
	    score_matrix.trunc_gap_length,
	    score_matrix.decay_gap_length);

	std::cout << std::string(buf) << std::endl;

	int n_input_seqs = 0;
	int max_input_seqs = 0;
	LPOSequence_T *seq=NULL,*lpo_out=NULL,*frame_seq=NULL,*dna_lpo=NULL,*lpo_in=NULL;
	LPOSequence_T **input_seqs=NULL;
	CALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
	char *print_matrix_letters=NULL,*fasta_out=NULL,*po_out=NULL,*clustal_out=NULL,*matrix_filename=NULL,
    *seq_filename=NULL,*frame_dna_filename=NULL,*po_filename=NULL,*po2_filename=NULL,
    *po_list_filename=NULL, *hbmin=NULL,*numeric_data=NULL,*numeric_data_name="Nmiscall",
    *dna_to_aa=NULL,*pair_score_file=NULL,*aafreq_file=NULL,*termval_file=NULL,
    *bold_seq_name=NULL,*subset_file=NULL,*subset2_file=NULL,*rm_subset_file=NULL,
    *rm_subset2_file=NULL;
    po_filename = NULL;
    float bundling_threshold=0.9;
    int exit_code=0,count_sequence_errors=0,please_print_snps=0,
    report_consensus_seqs=0,report_major_allele=0,use_aggressive_fusion=0;
    int show_allele_evidence=0,please_collapse_lines=0,keep_all_links=0;
    int remove_listed_seqs=0,remove_listed_seqs2=0,please_report_similarity;
    int do_global=1, do_progressive=0, do_preserve_sequence_order=0;
    int nseq=0,do_switch_case=dont_switch_case,do_analyze_bundles=0;
	//lpo_in = read_partial_order_file(po_filename, subset_file, remove_listed_seqs, keep_all_links, do_switch_case, &score_matrix);
	FILE *errfile=stderr,*logfile=NULL,*lpo_file_out=NULL,*po_list_file=NULL,*seq_ifile=NULL;


	do_switch_case=switch_case_to_lower;

	n_input_seqs = 0;
  	max_input_seqs = 10;
  	CALLOC (input_seqs, max_input_seqs, LPOSequence_T *);

	std::cout << "I did it\n";


	
	std::vector<std::string> input_strings;
	input_strings.push_back("ACATTAAATTAA\0\n");
	input_strings.push_back("GGGTTGGGTTGG\0\n");
	input_strings.push_back("GGGTTGCGTTGG\0\n");
	input_strings.push_back("AAATTAAATTAA\0\n");

	int isid = 0;
	do_switch_case = 0;
	char* test_name;
	char* test_title;
	char* new_seq;

	for (auto is : input_strings)
	{
		
    	std::string tns("test_name_"+std::to_string(isid));
    	std::string tts("test_title_"+std::to_string(isid));
		test_name = &tns[0];
		test_title = &tts[0];
		new_seq = &is[0];

		std::cout << "NAME:" << test_name << " TITLE: " << test_title << " SEQ: " << new_seq << "\n"; 
		
		create_seq(isid, &seq, test_name, test_title, new_seq, do_switch_case);
		isid += 1;
		nseq += 1;
	}



	for (int i=0; i<nseq; i++) {
		input_seqs[n_input_seqs++] = &(seq[i]);
		initialize_seqs_as_lpo(1,&(seq[i]),&score_matrix);
		if (n_input_seqs == max_input_seqs) {
			max_input_seqs *= 2;
			REALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
		}
	}

	/*
	seq_filename = "testifasta.fasta";

	if (seq_filename) {
    seq_ifile = fopen (seq_filename, "r");
    nseq = read_fasta (seq_ifile, &seq, do_switch_case, &comment);
    fclose (seq_ifile);

	std::cout << "I did it 2\n";

	//create_seq(nseq,seq,seq_name,seq_title,tmp_seq.p,do_switch_case
	*/

	/*
	seq_ifile = fopen(seq_filename, "r");

	nseq = read_fasta(seq_ifile, &seq, do_switch_case, &comment);

	fclose (seq_ifile);
	for (i=0; i<nseq; i++) {
		input_seqs[n_input_seqs++] = &(seq[i]);
		initialize_seqs_as_lpo(1,&(seq[i]),&score_matrix);
		if (n_input_seqs == max_input_seqs) {
			max_input_seqs *= 2;
			REALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
		}
	}

	*/

	std::cout << "I did it 2.5\n";
	
	lpo_out = buildup_progressive_lpo (isid, input_seqs, &score_matrix,
		use_aggressive_fusion, do_progressive, pair_score_file,
		matrix_scoring_function, do_global, do_preserve_sequence_order);

	std::cout << "I did it 3\n";

	generate_lpo_bundles(lpo_out,bundling_threshold);

	std::cout << "START\n";

	for (int i = 0; i < lpo_out->nsource_seq; i+=1)
	{
		std::cout << lpo_out->source_seq[i].name << "\n";
		std::cout << lpo_out->source_seq[i].title << "\n";
		std::cout << lpo_out->source_seq[i].length << "\n";
	    std::cout << lpo_out->source_seq[i].istart << "\n";
	    std::cout << lpo_out->source_seq[i].weight << "\n";
	    std::cout << lpo_out->source_seq[i].bundle_id << "\n";
	    //std::cout << lpo_out->source_seq[i].sequence << "\n";

	    std::cout << "---------------------------------------\n";
	}


	//std::cout << score_matrix.nsymbol << "\n";

	for (int i = 0; i < lpo_out->length; i+=1)
	{
		std::cout << score_matrix.symbol[lpo_out->letter[i].letter];
	}
	std::cout << "\n";


	char FASTA_GAP_CHARACTER = '-';
	int i,j,nring=0,iprint;
	int ibundle = -1;
	char **seq_pos=NULL, *p=NULL, *include_in_save=NULL;
	nring=xlate_lpo_to_al(lpo_out, score_matrix.nsymbol, score_matrix.symbol, ibundle, FASTA_GAP_CHARACTER, &seq_pos, &p, &include_in_save);


	std::cout << "START PRINTING SEQUENCES\n";
	std::cout << "******************************************\n";

	for (int i = 0; i < lpo_out->nsource_seq; i+=1)
	{
		std::cout << "SEQ NAME: " << lpo_out->source_seq[i].name << "\n";
		std::cout << "SEQ TITLE: " << lpo_out->source_seq[i].title << "\n";
		std::cout << "SEQ BUNDLE: " << lpo_out->source_seq[i].bundle_id << "\n";
		for (int j = 0; j < nring; j+=1)
		{
			std::cout << seq_pos[i][j];
		}
		std::cout << "\n";
		
	    //std::cout << lpo_out->source_seq[i].sequence << "\n";

	    std::cout << "******************************************\n";
	}

	
	//std::cout << type_name<decltype(typeid(lpo_out->length).name())>() << '\n';
	//std::cout << typeid(lpo_out->length).name() << "\n";
	//std::cout << lpo_out->length << "\n";
	//std::cout << lpo_out->title << "\n";
	//std::cout << lpo_out->sequence << "\n";
	//std::cout << lpo_out->letter[2].letter << "\n";

	/*
	std::cout << seq->letter[i].letter < score_matrix->nsymbol ? 
	    score_matrix->symbol[seq->letter[i].letter] 
	    : seq->letter[i].letter << "\n";
	*/


	std::cout << (*lpo_out).length << "\n";
	std::cout << (*lpo_out).title << "\n";
	std::cout << (*lpo_out).sequence << "\n";
	std::cout << (*lpo_out).letter->letter << "\n";

	std::cout << "END\n";

	//return;


	po_out = "poa_testi_output.lpo";
	clustal_out = "poa_testi_output.clustal";



	std::cout << "I did it 4\n";

	if (lpo_file_out=fopen(po_out, "w")) {
		write_lpo(lpo_file_out,lpo_out,&score_matrix);
		fclose(lpo_file_out);

		
	}

	if (seq_ifile=fopen(clustal_out,"w")) {

		export_clustal_seqal(seq_ifile,lpo_out,score_matrix.nsymbol,score_matrix.symbol);
		fclose(seq_ifile);

	}

	std::cout << "I did it 5\n";

	//std::function<LPOScore_T(int, int)>
	//std::function<int(int, int)>
	
}