#include "file_reader.hpp"

extern "C" 
{
 	#include "lpo.h"
	#include "msa_format.h"
	#include "align_score.h"
}

// Split the sequences into num_anchors+1 parts, construct POA separately for each part, and combine the split POAs into a single POA
// Arguments the same as for buildup_progressive_lpo except:
// - anchors is a two-dimensional num_anchors*nseq array where entry anchors[i][j] gives the position of the i:th anchor in the j:th sequence
// - num_anchors is the number of anchors
LPOSequence_T *split_align_combine(int nseq, LPOSequence_T **seqs, int **anchors, int num_anchors,
				   ResidueScoreMatrix_T *score_matrix,
				   int use_aggressive_fusion,
				   int do_progressive,
				   char score_file[], 
				   LPOScore_T (*scoring_function)
				   (int,int,LPOLetter_T [],LPOLetter_T [],
				    ResidueScoreMatrix_T *),
				   int use_global_alignment,
				   int preserve_sequence_order) {
  // Split the sequences according to anchors and do POA
  LPOSequence_T ***lpo_split;
  LPOSequence_T **lpo_out;
  int do_switch_case = 1;

  CALLOC(lpo_out, num_anchors+1, LPOSequence_T *);
  
  CALLOC(lpo_split, num_anchors+1, LPOSequence_T **);
  for(int i = 0; i < num_anchors+1; i++) {
    CALLOC(lpo_split[i], nseq, LPOSequence_T *);
    for(int j = 0; j < nseq; j++) {
      CALLOC(lpo_split[i][j], 1, LPOSequence_T);
      int end = i >= num_anchors ? seqs[j]->length : anchors[i][j];
      int start = i == 0 ? 0 : anchors[i-1][j];
      CALLOC(lpo_split[i][j]->sequence, end-start+1, char);
      lpo_split[i][j]->length = end-start;
      CALLOC(lpo_split[i][j]->title, SEQUENCE_NAME_MAX, char);
      strncpy(lpo_split[i][j]->title, seqs[j]->title, SEQUENCE_NAME_MAX);
      strncpy(lpo_split[i][j]->sequence, &seqs[j]->sequence[start], end-start);
      lpo_split[i][j]->sequence[end-start] = '\0';
      strncpy(lpo_split[i][j]->name, seqs[j]->name, SEQUENCE_NAME_MAX);
      lpo_split[i][j]->nsource_seq = 0;
      
      lpo_init(lpo_split[i][j]);
    }
    
    lpo_out[i] = buildup_progressive_lpo (nseq, lpo_split[i], score_matrix, 
					  use_aggressive_fusion, do_progressive, score_file, 
					  scoring_function, use_global_alignment, preserve_sequence_order);
  }

  // Combine the POA of split sequences

  // The length of the combined POA
  int newlen = 0;
  for(int i = 0; i < num_anchors+1; i++) {
    newlen += lpo_out[i]->length;
  }

  //std::cout << "New len: " << newlen << std::endl;

  // The combined POA
  LPOSequence_T *lpo_final;

  CALLOC(lpo_final, 1, LPOSequence_T);
  lpo_final->length = newlen;
  CALLOC(lpo_final->letter, newlen, LPOLetter_T);
  CALLOC(lpo_final->title, 10, char);
  strncpy(lpo_final->title, "FINAL", 10);
  CALLOC(lpo_final->sequence, newlen, char);
  strncpy(lpo_final->name, "FINAL", 10);
  lpo_final->nsource_seq = 0;
  
  // copy LPOLetters
  int offset = 0;          /* How much do we need to shift the split POA which we are currently copying */
  LPOLetterRef_T *posmap;  /* Mapping from old positions in split POA to new positions in combined POA */

  CALLOC(posmap, newlen, LPOLetterRef_T);

  for(int i = 0; i < num_anchors+1; i++) {
    // Set up the position mapping for current split POA
    for(int j = 0; j < lpo_out[i]->length; j++) {
      posmap[j] = offset+j;
    }
	
    for(int j = 0; j < lpo_out[i]->length; j++) {
      LPOLetter_T *fi_letter = &(lpo_final->letter[offset+j]);
      LPOLetter_T *old_letter = &(lpo_out[i]->letter[j]);

      // copy left
      LPOLetterLink_T *fi_next = &(fi_letter->left);
      LPOLetterLink_T *old_next = &(old_letter->left);
      while(old_next) {
	fi_next->ipos = old_next->ipos == -1 ? -1 : posmap[old_next->ipos];
	fi_next->more = NULL;
	
	old_next = old_next->ipos == -1 ? NULL : old_next->more;
	if (old_next) {
	  LPOLetterLink_T *tmp;
	  CALLOC(tmp, 1, LPOLetterLink_T);
	  fi_next->more = tmp;
	  fi_next = tmp;
	}
      }

      // copy right
      fi_next = &(fi_letter->right);
      old_next = &(old_letter->right);
      while(old_next) {
	fi_next->ipos = old_next->ipos == -1 ? -1 : posmap[old_next->ipos];
	fi_next->more = NULL;
	
	old_next = old_next->more;
	if (old_next) {
	  LPOLetterLink_T *tmp;
	  CALLOC(tmp, 1, LPOLetterLink_T);
	  fi_next->more = tmp;
	  fi_next = tmp;
	}
      }

      // copy source
      LPOLetterSource_T *fi_next2 = &(fi_letter->source);
      LPOLetterSource_T *old_next2 = &(old_letter->source);
      while(old_next2) {
	fi_next2->iseq = old_next2->iseq;
	if (old_next2->iseq >= 0) {
	  fi_next2->ipos = i == 0 ? old_next2->ipos : old_next2->ipos + anchors[i-1][old_next2->iseq];
	}
	fi_next2->more = NULL;
	
	old_next2 = old_next2->more;
	if (old_next2) {
	  LPOLetterSource_T *tmp;
	  CALLOC(tmp, 1, LPOLetterSource_T);
	  fi_next2->more = tmp;
	  fi_next2 = tmp;
	}
      }

      
      // copy align_ring
      fi_letter->align_ring = posmap[old_letter->align_ring];
      
      // copy ring_id
      fi_letter->ring_id = posmap[old_letter->ring_id];

      // copy score
      fi_letter->score = old_letter->score;
      
      // copy letter
      fi_letter->letter = old_letter->letter;

    }
    strncpy(&(lpo_final->sequence[offset]), lpo_out[i]->sequence, lpo_out[i]->length);
    
    offset += lpo_out[i]->length;
  }

  FREE(posmap);

  // Add links between the split LPOs
  offset = 0;
  for(int i = 0; i < num_anchors; i++) {
    for(int j = 0; j < lpo_out[i]->length; j++) {
      if (lpo_out[i]->letter[j].right.ipos == -1) {
	// offset+j is right end (i.e no further links to the right), find link in the next split POA
	int offset2 = offset + lpo_out[i]->length;  /* This is the offset of the next split POA */
	for(int k = 0; k < lpo_out[i+1]->length; k++) {
	  if (lpo_out[i+1]->letter[k].left.ipos == -1) {
	    // offset2+k is a left end (i.e. no further links to the left), see if it should be linked to offset+j (i.e. they share a source sequence)
	    LPOLetterSource_T *j_source=&(lpo_final->letter[j+offset].source);
	    LPOLetterSource_T *k_source=&(lpo_final->letter[k+offset2].source);
	    int ok = 0;
	    while(!ok && j_source) {
	      while(!ok && k_source) {
		if (j_source->iseq != -1 && k_source->iseq != -1 && j_source->iseq == k_source->iseq) {
		  // A shared source sequence found, establish the link
		  //std::cout << "Linking " << (j+offset) << " to " << (k+offset2) << std::endl;
										    
		  LPOLetterLink_T *j_link = &(lpo_final->letter[j+offset].right);
		  if (j_link->ipos == -1) {
		    j_link->ipos = k+offset2;
		    j_link->more = NULL;
		  } else {
		    while(j_link->more) j_link = j_link->more;
		    CALLOC(j_link->more, 1, LPOLetterLink_T);
		    j_link->more->ipos = k+offset2;
		    j_link->more->more = NULL;
		  }

		  
		  LPOLetterLink_T *k_link = &(lpo_final->letter[k+offset2].left);
		  if (k_link->ipos == -1) {
		    k_link->ipos = j+offset;
		    k_link->more = NULL;
		  } else {
		    while(k_link->more) k_link = k_link->more;
		    CALLOC(k_link->more, 1, LPOLetterLink_T);
		    k_link->more->ipos = j+offset;
		    k_link->more->more = NULL;
		  }
		  ok = 1;
		}
		k_source = k_source->more;
	      }
	      j_source = j_source->more;
	    }
	  }
	}
      }
    }
    offset += lpo_out[i]->length;
  }
  
  // Copy LPOSourceInfo
  for(int i = 0; i < nseq; i++) {
    save_lpo_source(lpo_final, seqs[i]->name, seqs[i]->title, seqs[i]->length, 1, NO_BUNDLE,0, NULL);
  }

  // Free memory
  for(int i = 0; i < num_anchors+1; i++) {
    free_lpo_sequence(lpo_out[i],0);
    for(int j = 0; j < nseq; j++) {
      free_lpo_sequence(lpo_split[i][j],1);
    }
    FREE(lpo_split[i]);
  }
  FREE(lpo_split);

  FREE(lpo_out);
      
  return lpo_final;
}

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


int main(int argc, char *argv[]) {
  // Sequences to align
  char *seqs[] = {
    "acgtgatagctagctagactagatcgacgagcagcatcgagggactacgactagcaggactatcatctaaggcatcacgacgactgcg",
    "acgtgatagctacgtagactagatcgacacgtataagcatcgagggactacgactagcaggactatcatctaacgtgcatcacgacgactacg",
    "acgtgatagcgttagctagactagaacactcgacgagcacacatcgagggactacgactagcagacgtatatatcatctaagcacgacgacgactaag",
    "acgtacgatagctagctagactagaccgacgagcagcatcgagggactacgaacagcaggactaacgatctaaggcatcacgtccgacgactacg"
  };

  int **anchors = new int *[3];
  for(int i = 0; i < 3; i++) {
    anchors[i] = new int[4];
  }

  anchors[0][0] = 6;
  anchors[0][1] = 6;
  anchors[0][2] = 6;
  anchors[0][3] = 8;

  anchors[1][0] = 40;
  anchors[1][1] = 43;
  anchors[1][2] = 47;
  anchors[1][3] = 42;

  anchors[2][0] = 79;
  anchors[2][1] = 84;
  anchors[2][2] = 89;
  anchors[2][3] = 86;

  //int do_switch_case=dont_switch_case;
  int do_switch_case=1;

  // Build score matrix
  ResidueScoreMatrix_T score_matrix;
  std::string score_matrix_path = "blosum80.mat";
  int path_len = score_matrix_path.length();
  char score_matrix_path_ca[path_len+1];
  strcpy(score_matrix_path_ca, score_matrix_path.c_str());
  int score_matrix_int = read_score_matrix(score_matrix_path_ca, &score_matrix);

  int seq_name_max = 50;
  int seq_title_max = 50;
  int seq_length_max = 100;

  char *seq_line;
  char *seq_title;

  int length = 0;
  int nseq = 4;
  int isid = 0;
	
  std::string name_str = "";
  std::string title_str = "";
  std::string line_str = "";
	
  // Use strings to create lpo
  LPOSequence_T **input_seqs = NULL;
	
  int n_input_seqs = 0;
  int max_input_seqs = 40;

  CALLOC (input_seqs, max_input_seqs, LPOSequence_T *);

  int di = 0;
  for (di=0; di<nseq; di++) {
    CALLOC(input_seqs[di], 1, LPOSequence_T);

    name_str = "tn" + std::to_string(di);
    title_str = "tt" + std::to_string(di);
    CALLOC(input_seqs[di]->title, seq_title_max, char);
    CALLOC(input_seqs[di]->sequence, seq_length_max, char);
    
    strcpy(input_seqs[di]->title, title_str.c_str());
    strcpy(input_seqs[di]->sequence, seqs[di]);
    strcpy(input_seqs[di]->name, name_str.c_str());
    input_seqs[di]->nsource_seq = 0;
    input_seqs[di]->length = strlen(seqs[di]);
      
    lpo_init(input_seqs[di]);
    
    if (n_input_seqs == max_input_seqs) {
      max_input_seqs *= 2;
      REALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
    }
    isid++;
  }

  
  int use_aggressive_fusion = 0;
  int do_progressive = 0;
  int do_global = 1;
  int do_preserve_sequence_order = 0;
  char * pair_score_file = NULL;

  LPOSequence_T *lpo_out = NULL;

  // Build progressive lpo
  //lpo_out = buildup_progressive_lpo (isid, input_seqs, &score_matrix,
  //		use_aggressive_fusion, do_progressive, pair_score_file,
  //		matrix_scoring_function, do_global, do_preserve_sequence_order);
  
  lpo_out = split_align_combine(isid, input_seqs, anchors, 3, &score_matrix, 
				use_aggressive_fusion, do_progressive, pair_score_file, 
				matrix_scoring_function, do_global, do_preserve_sequence_order);
  
  // Create lpo bundles
  double bundling_threshold = 0.8;
  generate_lpo_bundles_new(lpo_out, bundling_threshold);  


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

  for(int i = 0; i < 3; i++) {
    delete [] anchors[i];
  }
  delete [] anchors;
  
  FREE(p); 
  FREE(include_in_save);
  FREE(seq_pos);

  for (int di=0; di<nseq; di++) {
    free_lpo_sequence(input_seqs[di],1);
  }
  FREE(input_seqs);

  free_lpo_sequence(lpo_out,1);


  return 0;
}
