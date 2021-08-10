#include "fun_poa.hpp"
#include <iostream>
#include <vector>
#include <stdio.h>
#include <string>
#include <typeinfo>


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
  if (p_seq_pos)
    *p_seq_pos = seq_pos;
  return nring;
}


void LOMEXPOAGlue::initialize_me(std::string score_matrix_path)
{
	std::cout << "INIT START\n";
	int path_len = score_matrix_path.length();
	char score_matrix_path_ca[path_len+1];
	strcpy(score_matrix_path_ca, score_matrix_path.c_str());
	score_matrix_int = read_score_matrix(score_matrix_path_ca, &score_matrix);
	std::cout << "INIT END\n";
}


std::vector<std::string> LOMEXPOAGlue::poa_consensus_sequences(std::vector<std::string> input_sequences, float bundling_threshold)
{

	std::cout << "GLUE 1\n";
	std::vector<std::string> consensus_sequences;
	
	LPOSequence_T **input_seqs = NULL;
	LPOSequence_T *seq = NULL;
	LPOSequence_T *lpo_out = NULL;
	int n_input_seqs = 0;
	int max_input_seqs = 10;
	int nseq = 0;
	int isid = 0;
	int inseqnum = input_sequences.size();
	
	// Options ??
	int do_switch_case = 1; /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WHAT SHOULD I PUT HERE??
	int use_aggressive_fusion = 0;
	int do_progressive = 0;
	int do_global = 1;
	int do_preserve_sequence_order = 1;
	char * pair_score_file = NULL;



	CALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
			

	char * name;
	char * title;
	char * sequ;

	std::cout << "GLUE 2\n";

	// Read input sequences
	for (auto seq_str : input_sequences)
	{
    	std::string name_str("test_name_"+std::to_string(isid));
    	std::string title_str("test_title_"+std::to_string(isid));
		name = &name_str[0];
		title = &title_str[0];
		sequ = &seq_str[0];
		//std::cout << "NAME:" << test_name << " TITLE: " << test_title << " SEQ: " << new_seq << "\n"; 
		create_seq(isid, &seq, name, title, sequ, do_switch_case);
		isid += 1;
		nseq += 1;
	}

	std::cout << "GLUE 3\n";

	// Make LPO? out of the input sequences
	for (int i = 0; i < nseq; i+=1) {
		input_seqs[n_input_seqs++] = &(seq[i]);
		initialize_seqs_as_lpo(1,&(seq[i]),&score_matrix);
		if (n_input_seqs == max_input_seqs) {
			max_input_seqs *= 2;
			REALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
		}
	}


	std::cout << "GLUE 4\n";

	// Bundle consensus sequences
	lpo_out = buildup_progressive_lpo (isid, input_seqs, &score_matrix, 
		use_aggressive_fusion, do_progressive, pair_score_file, 
		matrix_scoring_function, do_global, do_preserve_sequence_order);
	

	std::cout << "GLUE 5\n";


	generate_lpo_bundles(lpo_out, bundling_threshold);


	std::cout << "GLUE 6\n";

	char gap_char = '-';
	int nring = 0;
	int ibundle = 1;
	// What are these?
	char **seq_pos = NULL;
	char *p = NULL;
	char *include_in_save = NULL;
	
	// xlate means ???
	nring = xlate_lpo_to_al(lpo_out, score_matrix.nsymbol, score_matrix.symbol, ibundle, gap_char, &seq_pos, &p, &include_in_save);


	std::cout << "GLUE 7\n";

	// Extract sequences
	std::cout << "****************** CONSENSUS SEQUENCES ************************\n";
	for (int i = inseqnum; i < lpo_out->nsource_seq; i += 1)
	{
		std::cout << "SEQ NAME: " << lpo_out->source_seq[i].name << "\n";
		std::cout << "SEQ TITLE: " << lpo_out->source_seq[i].title << "\n";
		std::cout << "SEQ BUNDLE: " << lpo_out->source_seq[i].bundle_id << "\n";
		for (int j = 0; j < nring; j+=1)
		{
			std::cout << seq_pos[i][j];
		}
		std::cout << "\n";
	    std::cout << "**************************************************************\n";
	}

	std::cout << "GLUE 8\n";

	// Free allocated memory
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
	FREE (input_seqs);
	if (nseq > 0){FREE (seq);}
	//FREE (lpo_out); /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DO I NEED THIS?

	// Return consensus sequences
	return consensus_sequences;
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

	int bbb = mini_test_from_align_score();

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