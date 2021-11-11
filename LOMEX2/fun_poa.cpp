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




//
//
//  New version that hopefully fixes some errors with bad inputs
//  
//  AND also even newer version that fixes how split graphs are combined
//
//  AND hopefully faster version WITH flipped for-for loop order 
//
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
  int *anchor_good_seqs;
  CALLOC(anchor_good_seqs, num_anchors+1, int);

  // Switches to consider all nodes from both combined graphs instead of the ones with no nodes on right/left
  bool DO_ALL_NODES_LEFT = true;
  bool DO_ALL_NODES_RIGHT = true;

  int do_switch_case = 1;

  CALLOC(lpo_out, num_anchors+1, LPOSequence_T *);
  
  CALLOC(lpo_split, num_anchors+1, LPOSequence_T **);

  int bad_anchor_blocks = 0;

  int nof_legit_anchors = num_anchors;

  for(int i = 0; i < num_anchors+1; i++) {

    int legit_anchor_index = i - bad_anchor_blocks;
    // Calculate how many non-empty blocks there are for this anchor
    int nof_legit_blocks = 0;
    for(int j = 0; j < nseq; j++) {
      int end = i >= num_anchors ? seqs[j]->length : anchors[i][j];
      int start = i == 0 ? 0 : anchors[i-1][j];
      if (end-start > 0){nof_legit_blocks+=1;}
    }

    anchor_good_seqs[i] = nof_legit_blocks;

    if (nof_legit_blocks == 0){
      bad_anchor_blocks += 1;
      REALLOC(lpo_out, num_anchors+1-bad_anchor_blocks, LPOSequence_T *);
      nof_legit_anchors -= 1;
      continue;
    }


    CALLOC(lpo_split[legit_anchor_index], nof_legit_blocks, LPOSequence_T *);

    int legit_block_seq_index = 0;

    for(int j = 0; j < nseq; j++) {

      
      int end = i >= num_anchors ? seqs[j]->length : anchors[i][j];
      int start = i == 0 ? 0 : anchors[i-1][j];

      if (end-start <= 0){continue;}

      CALLOC(lpo_split[legit_anchor_index][legit_block_seq_index], 1, LPOSequence_T);
      CALLOC(lpo_split[legit_anchor_index][legit_block_seq_index]->sequence, end-start+1, char);
      lpo_split[legit_anchor_index][legit_block_seq_index]->length = end-start;
      CALLOC(lpo_split[legit_anchor_index][legit_block_seq_index]->title, SEQUENCE_NAME_MAX, char);
      strncpy(lpo_split[legit_anchor_index][legit_block_seq_index]->title, seqs[j]->title, SEQUENCE_NAME_MAX);
      strncpy(lpo_split[legit_anchor_index][legit_block_seq_index]->sequence, &seqs[j]->sequence[start], end-start);
      lpo_split[legit_anchor_index][legit_block_seq_index]->sequence[end-start] = '\0';
      strncpy(lpo_split[legit_anchor_index][legit_block_seq_index]->name, seqs[j]->name, SEQUENCE_NAME_MAX);
      lpo_split[legit_anchor_index][legit_block_seq_index]->nsource_seq = 0;
      
      lpo_init(lpo_split[legit_anchor_index][legit_block_seq_index]);

      legit_block_seq_index += 1;
    }
    
    lpo_out[legit_anchor_index] = buildup_progressive_lpo (nof_legit_blocks, lpo_split[legit_anchor_index], score_matrix, 
            use_aggressive_fusion, do_progressive, score_file, 
            scoring_function, use_global_alignment, preserve_sequence_order);
  }

  // Combine the POA of split sequences

  // The length of the combined POA
  int newlen = 0;
  for(int i = 0; i < nof_legit_anchors+1; i++) {
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

  for(int i = 0; i < nof_legit_anchors+1; i++) {
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
  //int total_hits = 0;

  // Add links between the split LPOs
                                                                                              // HERE A LIST OF SOLVED SEQS

  offset = 0;


  for(int i = 0; i < nof_legit_anchors; i++) {
    // Use this to count how many sequences are linked
    //int solved_links[nseq] = {0};
    int nof_solved_links = 0;     
    //for(int j = 0; j < lpo_out[i]->length; j++) {                                                              // MODDED TO START FROM tHE END
    
    //for(int j = lpo_out[i]->length - 1; j >= 0; j--) {
    for(int k = 0; k < lpo_out[i+1]->length; k++) {
         
      if (nof_solved_links >= nseq){break;}

      // ONLY NODES THAT DO NOT HAVE RIGHT NODE CHECKED if switch not on!!!
      //if ((lpo_out[i]->letter[j].right.ipos == -1)  || (DO_ALL_NODES_LEFT)){
      if ((lpo_out[i+1]->letter[k].left.ipos == -1)  || (DO_ALL_NODES_RIGHT)){

        // offset+j is right end (i.e no further links to the right), find link in the next split POA
        int offset2 = offset + lpo_out[i]->length;  /* This is the offset of the next split POA */
        //for(int k = 0; k < lpo_out[i+1]->length; k++) {
        for(int j = lpo_out[i]->length - 1; j >= 0; j--) { 

          if (nof_solved_links >= nseq){break;}

          // ONLY NODES THAT DO NOT HAVE LEFT NODE CHECKED if switch not on!!!
          //if ((lpo_out[i+1]->letter[k].left.ipos == -1)  || (DO_ALL_NODES_RIGHT)){
          if ((lpo_out[i]->letter[j].right.ipos == -1)  || (DO_ALL_NODES_LEFT)){  
            // offset2+k is a left end (i.e. no further links to the left), see if it should be linked to offset+j (i.e. they share a source sequence)
            LPOLetterSource_T *j_source=&(lpo_final->letter[j+offset].source);
            LPOLetterSource_T *k_source=&(lpo_final->letter[k+offset2].source);
            int ok = 0;
            //while(!ok && j_source) {
            while(j_source) { 
              // If link exists, exit
              if (ok) {break;}
              // Added this line below to reset the right side 
              k_source=&(lpo_final->letter[k+offset2].source);
              //while(!ok && k_source) {
              while(k_source) {
                // If link exists, exit
                if (ok) {break;}

                //std::cout << "???????????????????????????????????????\n";
                //std::cout << "Left node id: " << j+offset << "\n";
                //std::cout << "Right node id: " << k+offset2 << "\n";
                //std::cout << "???????????????????????????????????????\n";
                // Check if the nodes should actually be connected
                //int lposi = lpo_out[i]->letter[j].

                /*
                int lposi = j_source->ipos;
                int rposi = k_source->ipos;
                int lseqi = j_source->iseq;
                int rseqi = k_source->iseq;
                bool is_a_hit = true;

                if ((lseqi == rseqi) || true){
                  std::cout << "----------------------------\n";
                  std::cout << "Left pos is " << lposi << "\n";
                  std::cout << "Right pos is " << rposi << "\n";
                  std::cout << "Left seq is " << lseqi << "\n";
                  std::cout << "Right seq is " << rseqi << "\n";
                  std::cout << "----------------------------\n";
                }
                if (rposi-lposi != 1){is_a_hit = false;}
                */

                /*
                if (sseqi != -1){ 
                  for (int ci = lposi +1; ci < rposi; ci++){
                    std::cout << lpo_final->source_seq[sseqi].sequence;
                    if (lpo_final->source_seq[sseqi].sequence[ci] == '-'){is_a_hit = false;}
                  }
                } else {
                  is_a_hit = false;
                }
                */
                                /*
                if (is_a_hit){
                  std::cout << "HIT\n";
                } else {
                  std::cout << "MISS\n";
                }*/

                // Modified this condition so that it also checks that the left side character and the right side character are consecutive in the sequence they come from
                //if (is_a_hit && j_source->iseq != -1 && k_source->iseq != -1 && j_source->iseq == k_source->iseq) {
                if (k_source->ipos - j_source->ipos == 1 && j_source->iseq != -1 && k_source->iseq != -1 && j_source->iseq == k_source->iseq) {

                  //total_hits += 1;                  
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
                  // Link found, set ok flag
                  ok = 1;
                }
                k_source = k_source->more;
              }
              j_source = j_source->more;
            }

            // If nodes linked, check which sequences were linked through this link
            if (ok){
              LPOLetterSource_T *j_source_4pairing = &(lpo_final->letter[j+offset].source);
              LPOLetterSource_T *k_source_4pairing = &(lpo_final->letter[k+offset2].source);
              while(j_source_4pairing) { 
                k_source_4pairing = &(lpo_final->letter[k+offset2].source);
                while(k_source_4pairing) {
                  if (k_source_4pairing->ipos - j_source_4pairing->ipos == 1 && j_source_4pairing->iseq != -1 && k_source_4pairing->iseq != -1 && j_source_4pairing->iseq == k_source_4pairing->iseq) {
                    //solved_links[j_source_4pairing->iseq] = 1;
                    nof_solved_links = nof_solved_links + 1;
                  }
                  k_source_4pairing = k_source_4pairing->more;
                }
                j_source_4pairing = j_source_4pairing->more;
              }
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

  for (int i = 0; i < nof_legit_anchors+1; i+=1){free_lpo_sequence(lpo_out[i],0);}
  int delid = 0;
  for (int i = 0; i < num_anchors+1; i+=1){
    if (anchor_good_seqs[i] == 0){continue;}
    for(int j = 0; j < anchor_good_seqs[i]; j+=1){free_lpo_sequence(lpo_split[delid][j],1);}
    FREE (lpo_split[delid]);
    delid += 1;
  }
  FREE(lpo_split);
  FREE(lpo_out);
  FREE (anchor_good_seqs);


  /*
  int delid = 0;
  for(int i = 0; i < num_anchors+1; i++) {
  //for(int i = 0; i < nof_legit_anchors+1; i++) {
    if (i < nof_legit_anchors+1){free_lpo_sequence(lpo_out[i],0);}
    for(int j = 0; j < anchor_good_seqs[i]; j++) {
      free_lpo_sequence(lpo_split[delid][j],1);
    }
    if (anchor_good_seqs[i] > 0){delid+=1;}
    FREE(lpo_split[i]);
  }
  FREE(lpo_split);

  FREE(lpo_out);

  FREE (anchor_good_seqs);
  */
  
  //std::cout << "Total hits: " << total_hits << "\n";
  return lpo_final;
}






//
//
//  New version that hopefully fixes some errors with bad inputs
//  
//  AND also even newer version that fixes how split graphs are combined
//
// Split the sequences into num_anchors+1 parts, construct POA separately for each part, and combine the split POAs into a single POA
// Arguments the same as for buildup_progressive_lpo except:
// - anchors is a two-dimensional num_anchors*nseq array where entry anchors[i][j] gives the position of the i:th anchor in the j:th sequence
// - num_anchors is the number of anchors
LPOSequence_T *split_align_combine_WORKING_probably(int nseq, LPOSequence_T **seqs, int **anchors, int num_anchors,
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
  int *anchor_good_seqs;
  CALLOC(anchor_good_seqs, num_anchors+1, int);

  // Switches to consider all nodes from both combined graphs instead of the ones with no nodes on right/left
  bool DO_ALL_NODES_LEFT = true;
  bool DO_ALL_NODES_RIGHT = true;

  int do_switch_case = 1;

  CALLOC(lpo_out, num_anchors+1, LPOSequence_T *);
  
  CALLOC(lpo_split, num_anchors+1, LPOSequence_T **);

  int bad_anchor_blocks = 0;

  int nof_legit_anchors = num_anchors;

  for(int i = 0; i < num_anchors+1; i++) {

    int legit_anchor_index = i - bad_anchor_blocks;
    // Calculate how many non-empty blocks there are for this anchor
    int nof_legit_blocks = 0;
    for(int j = 0; j < nseq; j++) {
      int end = i >= num_anchors ? seqs[j]->length : anchors[i][j];
      int start = i == 0 ? 0 : anchors[i-1][j];
      if (end-start > 0){nof_legit_blocks+=1;}
    }

    anchor_good_seqs[i] = nof_legit_blocks;

    if (nof_legit_blocks == 0){
      bad_anchor_blocks += 1;
      REALLOC(lpo_out, num_anchors+1-bad_anchor_blocks, LPOSequence_T *);
      nof_legit_anchors -= 1;
      continue;
    }


    CALLOC(lpo_split[legit_anchor_index], nof_legit_blocks, LPOSequence_T *);

    int legit_block_seq_index = 0;

    for(int j = 0; j < nseq; j++) {

      
      int end = i >= num_anchors ? seqs[j]->length : anchors[i][j];
      int start = i == 0 ? 0 : anchors[i-1][j];

      if (end-start <= 0){continue;}

      CALLOC(lpo_split[legit_anchor_index][legit_block_seq_index], 1, LPOSequence_T);
      CALLOC(lpo_split[legit_anchor_index][legit_block_seq_index]->sequence, end-start+1, char);
      lpo_split[legit_anchor_index][legit_block_seq_index]->length = end-start;
      CALLOC(lpo_split[legit_anchor_index][legit_block_seq_index]->title, SEQUENCE_NAME_MAX, char);
      strncpy(lpo_split[legit_anchor_index][legit_block_seq_index]->title, seqs[j]->title, SEQUENCE_NAME_MAX);
      strncpy(lpo_split[legit_anchor_index][legit_block_seq_index]->sequence, &seqs[j]->sequence[start], end-start);
      lpo_split[legit_anchor_index][legit_block_seq_index]->sequence[end-start] = '\0';
      strncpy(lpo_split[legit_anchor_index][legit_block_seq_index]->name, seqs[j]->name, SEQUENCE_NAME_MAX);
      lpo_split[legit_anchor_index][legit_block_seq_index]->nsource_seq = 0;
      
      lpo_init(lpo_split[legit_anchor_index][legit_block_seq_index]);

      legit_block_seq_index += 1;
    }
    
    lpo_out[legit_anchor_index] = buildup_progressive_lpo (nof_legit_blocks, lpo_split[legit_anchor_index], score_matrix, 
            use_aggressive_fusion, do_progressive, score_file, 
            scoring_function, use_global_alignment, preserve_sequence_order);
  }

  // Combine the POA of split sequences

  // The length of the combined POA
  int newlen = 0;
  for(int i = 0; i < nof_legit_anchors+1; i++) {
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

  for(int i = 0; i < nof_legit_anchors+1; i++) {
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
  //int total_hits = 0;

  // Add links between the split LPOs
  																																														// HERE A LIST OF SOLVED SEQS

  offset = 0;
  for(int i = 0; i < nof_legit_anchors; i++) {
    // Use this to count how many sequences are linked
    //int solved_links[nseq] = {0};
    int nof_solved_links = 0;     
    //for(int j = 0; j < lpo_out[i]->length; j++) {                                                                 // MODDED TO START FROM tHE END
    for(int j = lpo_out[i]->length - 1; j >= 0; j--) { 
      if (nof_solved_links >= nseq){break;}

      // ONLY NODES THAT DO NOT HAVE RIGHT NODE CHECKED if switch not on!!!
      if ((lpo_out[i]->letter[j].right.ipos == -1)  || (DO_ALL_NODES_LEFT)){

        // offset+j is right end (i.e no further links to the right), find link in the next split POA
        int offset2 = offset + lpo_out[i]->length;  /* This is the offset of the next split POA */
        for(int k = 0; k < lpo_out[i+1]->length; k++) {
          if (nof_solved_links >= nseq){break;}

          // ONLY NODES THAT DO NOT HAVE LEFT NODE CHECKED if switch not on!!!
          if ((lpo_out[i+1]->letter[k].left.ipos == -1)  || (DO_ALL_NODES_RIGHT)){
            // offset2+k is a left end (i.e. no further links to the left), see if it should be linked to offset+j (i.e. they share a source sequence)
            LPOLetterSource_T *j_source=&(lpo_final->letter[j+offset].source);
            LPOLetterSource_T *k_source=&(lpo_final->letter[k+offset2].source);
            int ok = 0;
            //while(!ok && j_source) {
            while(j_source) { 
              // If link exists, exit
              if (ok) {break;}
              // Added this line below to reset the right side 
              k_source=&(lpo_final->letter[k+offset2].source);
              //while(!ok && k_source) {
              while(k_source) {
                // If link exists, exit
                if (ok) {break;}

                //std::cout << "???????????????????????????????????????\n";
                //std::cout << "Left node id: " << j+offset << "\n";
                //std::cout << "Right node id: " << k+offset2 << "\n";
                //std::cout << "???????????????????????????????????????\n";
                // Check if the nodes should actually be connected
                //int lposi = lpo_out[i]->letter[j].

                /*
                int lposi = j_source->ipos;
                int rposi = k_source->ipos;
                int lseqi = j_source->iseq;
                int rseqi = k_source->iseq;
                bool is_a_hit = true;

                if ((lseqi == rseqi) || true){
                  std::cout << "----------------------------\n";
                  std::cout << "Left pos is " << lposi << "\n";
                  std::cout << "Right pos is " << rposi << "\n";
                  std::cout << "Left seq is " << lseqi << "\n";
                  std::cout << "Right seq is " << rseqi << "\n";
                  std::cout << "----------------------------\n";
                }
                if (rposi-lposi != 1){is_a_hit = false;}
                */

                /*
                if (sseqi != -1){ 
                  for (int ci = lposi +1; ci < rposi; ci++){
                    std::cout << lpo_final->source_seq[sseqi].sequence;
                    if (lpo_final->source_seq[sseqi].sequence[ci] == '-'){is_a_hit = false;}
                  }
                } else {
                  is_a_hit = false;
                }
                */
                                /*
                if (is_a_hit){
                  std::cout << "HIT\n";
                } else {
                  std::cout << "MISS\n";
                }*/

                // Modified this condition so that it also checks that the left side character and the right side character are consecutive in the sequence they come from
                //if (is_a_hit && j_source->iseq != -1 && k_source->iseq != -1 && j_source->iseq == k_source->iseq) {
                if (k_source->ipos - j_source->ipos == 1 && j_source->iseq != -1 && k_source->iseq != -1 && j_source->iseq == k_source->iseq) {

                  //total_hits += 1;                  
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
                  // Link found, set ok flag
                  ok = 1;
                }
                k_source = k_source->more;
              }
              j_source = j_source->more;
            }

            // If nodes linked, check which sequences were linked through this link
            if (ok){
              LPOLetterSource_T *j_source_4pairing = &(lpo_final->letter[j+offset].source);
              LPOLetterSource_T *k_source_4pairing = &(lpo_final->letter[k+offset2].source);
              while(j_source_4pairing) { 
                k_source_4pairing = &(lpo_final->letter[k+offset2].source);
                while(k_source_4pairing) {
                  if (k_source_4pairing->ipos - j_source_4pairing->ipos == 1 && j_source_4pairing->iseq != -1 && k_source_4pairing->iseq != -1 && j_source_4pairing->iseq == k_source_4pairing->iseq) {
                    //solved_links[j_source_4pairing->iseq] = 1;
                    nof_solved_links = nof_solved_links + 1;
                  }
                }
              }
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

  for (int i = 0; i < nof_legit_anchors+1; i+=1){free_lpo_sequence(lpo_out[i],0);}
  int delid = 0;
  for (int i = 0; i < num_anchors+1; i+=1){
    if (anchor_good_seqs[i] == 0){continue;}
    for(int j = 0; j < anchor_good_seqs[i]; j+=1){free_lpo_sequence(lpo_split[delid][j],1);}
    FREE (lpo_split[delid]);
    delid += 1;
  }
  FREE(lpo_split);
  FREE(lpo_out);
  FREE (anchor_good_seqs);


  /*
  int delid = 0;
  for(int i = 0; i < num_anchors+1; i++) {
  //for(int i = 0; i < nof_legit_anchors+1; i++) {
    if (i < nof_legit_anchors+1){free_lpo_sequence(lpo_out[i],0);}
    for(int j = 0; j < anchor_good_seqs[i]; j++) {
      free_lpo_sequence(lpo_split[delid][j],1);
    }
    if (anchor_good_seqs[i] > 0){delid+=1;}
    FREE(lpo_split[i]);
  }
  FREE(lpo_split);

  FREE(lpo_out);

  FREE (anchor_good_seqs);
  */
  
  //std::cout << "Total hits: " << total_hits << "\n";
  return lpo_final;
}



//
//
//  New version that hopefully fixes some errors with bad inputs (has problems with how split graphs are combined)
//
//
// Split the sequences into num_anchors+1 parts, construct POA separately for each part, and combine the split POAs into a single POA
// Arguments the same as for buildup_progressive_lpo except:
// - anchors is a two-dimensional num_anchors*nseq array where entry anchors[i][j] gives the position of the i:th anchor in the j:th sequence
// - num_anchors is the number of anchors
LPOSequence_T *split_align_combine_oldfix(int nseq, LPOSequence_T **seqs, int **anchors, int num_anchors,
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
  int *anchor_good_seqs;
  CALLOC(anchor_good_seqs, num_anchors+1, int);



  int do_switch_case = 1;

  CALLOC(lpo_out, num_anchors+1, LPOSequence_T *);
  
  CALLOC(lpo_split, num_anchors+1, LPOSequence_T **);

  int bad_anchor_blocks = 0;

  int nof_legit_anchors = num_anchors;

  for(int i = 0; i < num_anchors+1; i++) {

    int legit_anchor_index = i - bad_anchor_blocks;
    // Calculate how many non-empty blocks there are for this anchor
    int nof_legit_blocks = 0;
    for(int j = 0; j < nseq; j++) {
      int end = i >= num_anchors ? seqs[j]->length : anchors[i][j];
      int start = i == 0 ? 0 : anchors[i-1][j];
      if (end-start > 0){nof_legit_blocks+=1;}
    }

    anchor_good_seqs[i] = nof_legit_blocks;

    if (nof_legit_blocks == 0){
      bad_anchor_blocks += 1;
      REALLOC(lpo_out, num_anchors+1-bad_anchor_blocks, LPOSequence_T *);
      nof_legit_anchors -= 1;
      continue;
    }


    CALLOC(lpo_split[legit_anchor_index], nof_legit_blocks, LPOSequence_T *);

    int legit_block_seq_index = 0;

    for(int j = 0; j < nseq; j++) {

      
      int end = i >= num_anchors ? seqs[j]->length : anchors[i][j];
      int start = i == 0 ? 0 : anchors[i-1][j];

      if (end-start <= 0){continue;}

      CALLOC(lpo_split[legit_anchor_index][legit_block_seq_index], 1, LPOSequence_T);
      CALLOC(lpo_split[legit_anchor_index][legit_block_seq_index]->sequence, end-start+1, char);
      lpo_split[legit_anchor_index][legit_block_seq_index]->length = end-start;
      CALLOC(lpo_split[legit_anchor_index][legit_block_seq_index]->title, SEQUENCE_NAME_MAX, char);
      strncpy(lpo_split[legit_anchor_index][legit_block_seq_index]->title, seqs[j]->title, SEQUENCE_NAME_MAX);
      strncpy(lpo_split[legit_anchor_index][legit_block_seq_index]->sequence, &seqs[j]->sequence[start], end-start);
      lpo_split[legit_anchor_index][legit_block_seq_index]->sequence[end-start] = '\0';
      strncpy(lpo_split[legit_anchor_index][legit_block_seq_index]->name, seqs[j]->name, SEQUENCE_NAME_MAX);
      lpo_split[legit_anchor_index][legit_block_seq_index]->nsource_seq = 0;
      
      lpo_init(lpo_split[legit_anchor_index][legit_block_seq_index]);

      legit_block_seq_index += 1;
    }
    
    lpo_out[legit_anchor_index] = buildup_progressive_lpo (nof_legit_blocks, lpo_split[legit_anchor_index], score_matrix, 
					  use_aggressive_fusion, do_progressive, score_file, 
					  scoring_function, use_global_alignment, preserve_sequence_order);
  }

  // Combine the POA of split sequences

  // The length of the combined POA
  int newlen = 0;
  for(int i = 0; i < nof_legit_anchors+1; i++) {
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

  for(int i = 0; i < nof_legit_anchors+1; i++) {
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
  for(int i = 0; i < nof_legit_anchors; i++) {
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

  for (int i = 0; i < nof_legit_anchors+1; i+=1){free_lpo_sequence(lpo_out[i],0);}
  int delid = 0;
  for (int i = 0; i < num_anchors+1; i+=1){
    if (anchor_good_seqs[i] == 0){continue;}
    for(int j = 0; j < anchor_good_seqs[i]; j+=1){free_lpo_sequence(lpo_split[delid][j],1);}
    FREE (lpo_split[delid]);
    delid += 1;
  }
  FREE(lpo_split);
  FREE(lpo_out);
  FREE (anchor_good_seqs);


  /*
  int delid = 0;
  for(int i = 0; i < num_anchors+1; i++) {
  //for(int i = 0; i < nof_legit_anchors+1; i++) {
    if (i < nof_legit_anchors+1){free_lpo_sequence(lpo_out[i],0);}
    for(int j = 0; j < anchor_good_seqs[i]; j++) {
      free_lpo_sequence(lpo_split[delid][j],1);
    }
    if (anchor_good_seqs[i] > 0){delid+=1;}
    FREE(lpo_split[i]);
  }
  FREE(lpo_split);

  FREE(lpo_out);

  FREE (anchor_good_seqs);
  */
      
  return lpo_final;
}


//
//
//  Old version which can cause errors with degenerate input
//
//
// Split the sequences into num_anchors+1 parts, construct POA separately for each part, and combine the split POAs into a single POA
// Arguments the same as for buildup_progressive_lpo except:
// - anchors is a two-dimensional num_anchors*nseq array where entry anchors[i][j] gives the position of the i:th anchor in the j:th sequence
// - num_anchors is the number of anchors
LPOSequence_T *split_align_combine_ORIGINAL(int nseq, LPOSequence_T **seqs, int **anchors, int num_anchors,
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
  CALLOC(lpo_split, num_anchors+1, LPOSequence_T **);
  
  LPOSequence_T **lpo_out;
  CALLOC(lpo_out, num_anchors+1, LPOSequence_T *);

  int do_switch_case = 1;

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

      //seq->letter[seq->length -1].right.ipos= INVALID_LETTER_POSITION;
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



// THE MAIN FUNCTION THAT IS USED TO CONSTRUCT CONSENSUS SEQUENCES
std::tuple<int, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, uint64_t> construct_consensus_seqs(int split_align, int aggro, int prog, int glob, int nof_anchors, int **anchors, int minkmerid, std::vector<std::string>& occurrence_sequences, std::vector<std::string>& consensus_sequences, 
																																																								ResidueScoreMatrix_T *score_matrix, int do_switch_case, float bundling_threshold, float support_threshold, 
																																																								int min_bundle_size_threshold, int min_char_support_threshold)
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
	int seq_length_max = 300;

	char seq_name[seq_name_max]="";
	char seq_line[seq_length_max]="";
	char seq_title[seq_title_max]="";

	stringptr tmp_seq=STRINGPTR_EMPTY_INIT;

	int length = 0;
	//int nseq = occurrence_sequences.size();
	int nseq = 0;
	int isid = 0;
	int inseqnum = occurrence_sequences.size();
	
	std::string name_str = "";
	std::string title_str = "";
	std::string line_str = "";
	

	poa_p1_start_time = std::chrono::high_resolution_clock::now();

	

// REPLACED BY ANOTHER VERSION

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

	


	// Use strings to create lpo
	
	/*
	LPOSequence_T **input_seqs = NULL;
	
	int n_input_seqs = 0;
	int max_input_seqs = 40;  //                                                       WAS 10!!!!!!!!!!!!!!

	CALLOC (input_seqs, max_input_seqs, LPOSequence_T *);

	int di = 0;
  for (di=0; di<nseq; di++) {
    CALLOC(input_seqs[di], 1, LPOSequence_T);

    name_str = "tn" + std::to_string(di);
    title_str = "tt" + std::to_string(di);
    CALLOC(input_seqs[di]->title, seq_title_max, char);
    CALLOC(input_seqs[di]->sequence, seq_length_max, char);

    strcpy(input_seqs[di]->title, title_str.c_str());
    //strcpy(input_seqs[di]->sequence, seqs[di]);
    strcpy(input_seqs[di]->sequence, occurrence_sequences[di].c_str());
    strcpy(input_seqs[di]->name, name_str.c_str());
    input_seqs[di]->nsource_seq = 0;
    //input_seqs[di]->length = strlen(seqs[di]);
    input_seqs[di]->length = strlen(occurrence_sequences[di].c_str());

    lpo_init(input_seqs[di]);

    if (n_input_seqs == max_input_seqs) {
      max_input_seqs *= 2;
      REALLOC (input_seqs, max_input_seqs, LPOSequence_T *);
    }
    isid++;
  }
	*/



	poa_p2_end_time = std::chrono::high_resolution_clock::now();


	// Bottleneck
	//====================================================================================
	poa_p3_start_time = std::chrono::high_resolution_clock::now();

	//int use_aggressive_fusion = 0; // CHANGED
	int use_aggressive_fusion = aggro;
	int do_progressive = prog;
	int do_global = glob; 																																																											// DO GLOBAL WAS PREVIOUSLY 0 !!!!!!!!!!!!!!!!!!!!!!
	//int do_global = 0; 
	int do_preserve_sequence_order = 0;
	char * pair_score_file = NULL;
	//float bundling_threshold = 0.4;

	LPOSequence_T *lpo_out = NULL;


	// Build progressive lpo

	if (split_align == 1){
		lpo_out = split_align_combine(isid, input_seqs, anchors, nof_anchors, score_matrix, 
				use_aggressive_fusion, do_progressive, pair_score_file, 
				matrix_scoring_function, do_global, do_preserve_sequence_order);
	} else {
		lpo_out = buildup_progressive_lpo (isid, input_seqs, score_matrix, 
				use_aggressive_fusion, do_progressive, pair_score_file, 
				matrix_scoring_function, do_global, do_preserve_sequence_order);
	}



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


	/*
  // THIS BLOCK PRINTS THE LPO GRAPH
	std::cout << "===============================================================================\n";
	std::cout << "HERE IS THE LPO GRAPH\n";
	write_lpo_modified(minkmerid, lpo_out, score_matrix);
	std::cout << "===============================================================================\n";
	//*/
	

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
		//*
		if (i < inseqnum){
			if (lpo_out->source_seq[i].bundle_id >= 0){
				std::cout << lpo_out->source_seq[i].bundle_id << "  " << conseq << "\n";
			} else {
				std::cout << lpo_out->source_seq[i].bundle_id << " " << conseq << "\n";	
			}
		} 
		//*/

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
			int supported_start = 0;
			for (int csi = 0; csi < nring; csi+=1)
			{
				if (character_supports.at(csi) == 1){supported_start = csi;break;}
			}

			// Find supported end
			int supported_end = 0;
			for (int csi = nring-1; csi >= 0; csi-=1)
			{
				if (character_supports.at(csi) == 1){supported_end = csi;break;}
			}
			
			// No random unsupported characters in the middle
			int from_support_to_no_support = 0;
			int from_no_support_to_support = 0;
			int illegal_dip_in_no_support = 0;
			int last_non_deletion = character_supports.at(0);
			for (int csi = 1; csi < nring; csi+=1)
			{
				if ((character_supports.at(csi) == -1) && (last_non_deletion == 1)){from_support_to_no_support+=1;}
				if ((character_supports.at(csi) == 1) && (last_non_deletion == -1)){
					from_no_support_to_support+=1;
					if (from_support_to_no_support > 0){illegal_dip_in_no_support += 1;}
				}
				if (character_supports.at(csi) == -1){last_non_deletion = -1;}
				if (character_supports.at(csi) == 1){last_non_deletion = 1;}
			}

			// Determine if enough support
			if (from_support_to_no_support > 1){enough_support = false;}
			if (from_no_support_to_support > 1){enough_support = false;}
			if (illegal_dip_in_no_support > 0){enough_support = false;}

			// Print accordingly
			if (enough_support && is_real_bundle && is_supported_percent && is_supported_hard){

				//consensus_sequences.push_back(conseq);
				consensus_sequences.push_back(conseq.substr(supported_start, supported_end-supported_start+1));
				std::cout << lpo_out->source_seq[i].bundle_id << "  " << conseq << " ** ("<< supported_start << "," << supported_end << ")\n";
			} else {
				std::cout << lpo_out->source_seq[i].bundle_id << "  " << conseq << " * \n";
				//continue;

				// SUS continue?



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

  // 																																								                                                        THIS ONLY USED WHEN SPLIT-ALIGN-METHOD IN USE
  if (split_align == 1){
  	free_lpo_sequence(lpo_out,1);
	}

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
