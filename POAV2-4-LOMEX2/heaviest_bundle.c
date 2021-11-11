

#include "default.h"
#include "poa.h"
#include "seq_util.h"
#include "lpo.h"


/* My modded version that should work */
/* Best path is chosen based on the already used path */
LPOLetterRef_T *heaviest_bundle_MY_NEWEST_MOD(int len,LPOLetter_T seq[],
        int nsource_seq,LPOSourceInfo_T source_seq[],
        int *p_best_len)
{
  int i,j,best_right,iright,ibest= -1,best_len=0;
  LPOLetterRef_T *best_path=NULL,*path=NULL;
  LPOLetterLink_T *right;
  LPOLetterSource_T *source;
  LPOScore_T *score=NULL,best_score= -999999,right_score;
  int *contains_pos=NULL,my_overlap,right_overlap;

  CALLOC(path,len,LPOLetterRef_T); /* GET MEMORY FOR DYNAMIC PROGRAMMING */
  CALLOC(score,len,LPOScore_T);
  CALLOC(contains_pos,nsource_seq,int);


  int initial_seq_weight = 1;
  int path_multiplier = 0;

  // ALLOCATE MEMORY FOR EACH INPUT SEQUENCE TO HOLD AN INTEGER WHICH TELLS HOW MANY TIMES IT HAS BEEN USED IN THE CURRENT CONSENSUS PATH
  
  // Initialize array for storing sequence counts in path stuff
  int **seq_in_path = (int**)malloc(len*sizeof(int *));
  for (int oko = 0; oko < len; oko++)
    seq_in_path[oko] = (int *)malloc(nsource_seq*sizeof(int));


  int *seq_appearing_in_path = NULL;
  CALLOC(seq_appearing_in_path, nsource_seq, int);

   int *seq_appearing_in_path_in_loop = NULL;
  CALLOC(seq_appearing_in_path_in_loop, nsource_seq, int);

  //printf("len is: %d \n", len);
  //printf("nsource_seq is: %d \n", nsource_seq);

  // INITIALIZE ALL TO ONE (ZERO WON'T WORK MAYBE?)
  LOOP (i,nsource_seq) /* RESET TO NOT MATCH ANY POSITIONS */
    contains_pos[i]= INVALID_LETTER_POSITION;

  // Initialize all to zero (ZERO WON'T WORK MAYBE?)
  for (int uro = 0; uro < len; uro++)
    for (int omnath = 0; omnath < nsource_seq; omnath++)
      seq_in_path[uro][omnath] = 0;

  LOOP (i,nsource_seq)
    seq_appearing_in_path[i] = initial_seq_weight;

  LOOP (i,nsource_seq)
    seq_appearing_in_path_in_loop[i] = initial_seq_weight;


  // For every node
  LOOPB (i,len) { /* FIND HEAVIEST PATH BY DYNAMIC PROGRAMMING */
    //printf("This is node %d \n", i);
    source= &seq[i].source;  /* MARK SEQUENCES CONTAINING THIS POSITION */
    memset(contains_pos,0,nsource_seq*sizeof(int)); /* ERASE ARRAY */

    // Empty the count for this node
    memset(seq_appearing_in_path,0,nsource_seq*sizeof(int)); /* ERASE ARRAY */


    // seq is a letter or in other words a node
    // seq[i] is the i:th node in the whole graph
    // seq[i].source is LPOLetterSource which is a linked list of sequence+position pairs which correspond to the node

    // source_seq is LPOSourceInfo list which has information about all the source/input sequences
    // source->iseq is source sequence id that corresponds to the current node

    do
      if (source_seq[source->iseq].weight>0){ /* EXCLUDE ZERO WEIGHT SEQS */
        contains_pos[source->iseq]=source->ipos+1; /* right MUST BE ADJACENT*/
        //printf("source->iseq is: %d \n", source->iseq);
        //printf("and contains pos for it is set to: %d \n", source->ipos+1);
      }
    while (source=source->more); /* KEEP COUNTING TILL NO more */

    // Now the whole source thing is exhausted i.e. we have looked at all sequences corresponding to the current node

    // Next

    right_score=right_overlap=0;  /*DEFAULT MOVE: NOTHING TO THE RIGHT*/
    best_right= INVALID_LETTER_POSITION;

    // Find best right?

    // Check all to the right (all valid rights?)
    for (right= &seq[i].right;right && right->ipos>=0;right=right->more) {
      my_overlap=0; /* OVERLAP CALCULATION */
      source= &seq[right->ipos].source;/*COUNT SEQS SHARED IN i AND right*/

      memset(seq_appearing_in_path_in_loop,0,nsource_seq*sizeof(int)); /* ERASE ARRAY */

      //printf("right->ipos i.e. nodes to the right: %d \n", right->ipos);

      // right is seq[i] node's linked node to the right 
      // right->ipos is probably an integer that is the id of the linked node
      // seq[right->ipos].source is LPOLetterSource which is a linked list of sequence+position pairs which correspond to the node
      // and now the source is the first sequence associated with the current linked node to the right

      
      do /* BIAS OVERLAP CALCULATION BY SEQUENCE WEIGHTING */
        if (contains_pos[source->iseq]==source->ipos){ /* YES, ADJACENT! */
          // My multiplier
          //printf("WE ARE ADJACENT\n");
          path_multiplier = seq_in_path[right->ipos][source->iseq];

          if (path_multiplier <= 0)
            path_multiplier = 1;

          // Adjust using the multiplier
          my_overlap += (source_seq[source->iseq].weight * path_multiplier);
          seq_appearing_in_path_in_loop[source->iseq] = 1;
        }
      while (source=source->more); /* KEEP COUNTING TILL NO more */


      //printf("DEBUG1\n");      

      /* FIND BEST RIGHT MOVE: BEST OVERLAP */
      if (my_overlap>right_overlap || (my_overlap==right_overlap && score[right->ipos]>right_score)) {
        right_overlap=my_overlap;
        right_score=score[right->ipos];
        best_right=right->ipos;
        // Update current best path sequences
        for (int gary = 0; gary < nsource_seq; gary++)
          seq_appearing_in_path[gary] = seq_appearing_in_path_in_loop[gary];
      }
    }

    //printf("DEBUG2\n");

    // Update counts for sequences in paths
    if (best_right >= 0){
    for (int hapatra = 0; hapatra < nsource_seq; hapatra++){
      if (seq_appearing_in_path[hapatra] == 1){
        seq_in_path[i][hapatra] = seq_in_path[best_right][hapatra] + 1;
      } else {
        seq_in_path[i][hapatra] = seq_in_path[best_right][hapatra];
      }
    }
    }

    //printf("DEBUG3\n");

    path[i]=best_right; /* SAVE THE BEST PATH FOUND */
    score[i]=right_score+right_overlap; /* SAVE THE SCORE */
    if (score[i]>best_score) { /* RECORD BEST SCORE IN WHOLE LPO */
      ibest=i;
      best_score=score[i];
    }
  }

  /*
  for (int rr = 0; rr < len; rr++){
    for (int cc = 0; cc < nsource_seq; cc++){
      printf("%d ", seq_in_path[rr][cc]);
    }
    printf("\n");
  }
  */

  CALLOC(best_path,len,LPOLetterRef_T); /* MEMORY FOR STORING BEST PATH */
  for (;ibest>=0;ibest=path[ibest])  /* BACK TRACK THE BEST PATH */
    best_path[best_len++]=ibest;

  FREE(path); /* DUMP SCRATCH MEMORY */
  FREE(score);
  FREE(contains_pos);

  for (int oko = 0; oko < len; oko++)
    FREE(seq_in_path[oko]);

  FREE(seq_in_path);
  FREE(seq_appearing_in_path);
  FREE(seq_appearing_in_path_in_loop);

  if (p_best_len) /* RETURN best_path AND ITS LENGTH */
    *p_best_len = best_len;
  return best_path;
}


/** finds the heaviest traversal of the LPO seq[], using dynamic programming;
  at each node the heaviest link is chosen to buildup traversals; finally,
  the traversal with the heaviest overall link weight is returned as an
  array of position indices.  The length of the array is stored in 
  *p_best_len*/

/* THE ORIGINAL HEAVIEST BUNDLE ALGORITHM */
LPOLetterRef_T *heaviest_bundle(int len,LPOLetter_T seq[],
				int nsource_seq,LPOSourceInfo_T source_seq[],
				int *p_best_len)
{
  int i,j,best_right,iright,ibest= -1,best_len=0;
  LPOLetterRef_T *best_path=NULL,*path=NULL;
  LPOLetterLink_T *right;
  LPOLetterSource_T *source;
  LPOScore_T *score=NULL,best_score= -999999,right_score;
  int *contains_pos=NULL,my_overlap,right_overlap;

  CALLOC(path,len,LPOLetterRef_T); /* GET MEMORY FOR DYNAMIC PROGRAMMING */
  CALLOC(score,len,LPOScore_T);
  CALLOC(contains_pos,nsource_seq,int);
  LOOP (i,nsource_seq) /* RESET TO NOT MATCH ANY POSITIONS */
    contains_pos[i]= INVALID_LETTER_POSITION;

  LOOPB (i,len) { /* FIND HEAVIEST PATH BY DYNAMIC PROGRAMMING */
    source= &seq[i].source;  /* MARK SEQUENCES CONTAINING THIS POSITION */
    memset(contains_pos,0,nsource_seq*sizeof(int)); /* ERASE ARRAY */
    do
      if (source_seq[source->iseq].weight>0) /* EXCLUDE ZERO WEIGHT SEQS */
	contains_pos[source->iseq]=source->ipos+1; /* right MUST BE ADJACENT*/
    while (source=source->more); /* KEEP COUNTING TILL NO more */

    right_score=right_overlap=0;  /*DEFAULT MOVE: NOTHING TO THE RIGHT*/
    best_right= INVALID_LETTER_POSITION;
    for (right= &seq[i].right;right && right->ipos>=0;right=right->more) {
      my_overlap=0; /* OVERLAP CALCULATION */
      source= &seq[right->ipos].source;/*COUNT SEQS SHARED IN i AND right*/
      do /* BIAS OVERLAP CALCULATION BY SEQUENCE WEIGHTING */
	if (contains_pos[source->iseq]==source->ipos) /* YES, ADJACENT! */
	  my_overlap += source_seq[source->iseq].weight;
      while (source=source->more); /* KEEP COUNTING TILL NO more */
      
      if (my_overlap>right_overlap /* FIND BEST RIGHT MOVE: BEST OVERLAP */
	  || (my_overlap==right_overlap && score[right->ipos]>right_score)) {
	right_overlap=my_overlap;
	right_score=score[right->ipos];
	best_right=right->ipos;
      }
    }

    path[i]=best_right; /* SAVE THE BEST PATH FOUND */
    score[i]=right_score+right_overlap; /* SAVE THE SCORE */
    if (score[i]>best_score) { /* RECORD BEST SCORE IN WHOLE LPO */
      ibest=i;
      best_score=score[i];
    }
  }

  CALLOC(best_path,len,LPOLetterRef_T); /* MEMORY FOR STORING BEST PATH */
  for (;ibest>=0;ibest=path[ibest])  /* BACK TRACK THE BEST PATH */
    best_path[best_len++]=ibest;

  FREE(path); /* DUMP SCRATCH MEMORY */
  FREE(score);
  FREE(contains_pos);

  if (p_best_len) /* RETURN best_path AND ITS LENGTH */
    *p_best_len = best_len;
  return best_path;
}




int assign_sequence_bundle_id(int path_length,LPOLetterRef_T path[],
			      LPOSequence_T *seq,int bundle_id,
			      float minimum_fraction)
{
  int i,*bundle_count=NULL,nseq_in_bundle=0;
  LPOLetterSource_T *source;
  
  CALLOC(bundle_count,seq->nsource_seq,int);
  LOOP (i,path_length) /* COUNT #POSITIONS OF EACH SEQ ARE IN path */
    for (source= &seq->letter[path[i]].source;source;source=source->more)
      bundle_count[source->iseq]++;

  LOOP (i,seq->nsource_seq) {/* FOR EACH SEQ OVER THRESHOLD, ASSIGN bundle_id*/
    /* printf("bundle %d:\t%s\t%d/%d %d\n",bundle_id,seq->source_seq[i].name, */
      /* bundle_count[i],seq->source_seq[i].length,seq->source_seq[i].weight); */
    if (seq->source_seq[i].bundle_id<0 /* NOT YET BUNDLED */
	&& seq->source_seq[i].length*minimum_fraction <= bundle_count[i]) {
      /* printf("   +++++++++++++++++"); */
      seq->source_seq[i].bundle_id = bundle_id; /* ASSIGN TO THIS BUNDLE */
      seq->source_seq[i].weight = 0; /* REMOVE FROM FUTURE heaviest_bundle */
      nseq_in_bundle++;
    }
    /* printf("\n"); */
  }

  FREE(bundle_count);
  return nseq_in_bundle; /* RETURN COUNT OF SEQUENCES IN BUNDLE */
}


int reweigh_sequences_in_bundle(int path_length,LPOLetterRef_T path[],
				LPOSequence_T *seq,int weight_multiplier,
				float minimum_fraction, int *reweighted)
{
  int i,*bundle_count=NULL,nseq_in_bundle=0;
  LPOLetterSource_T *source;
  
  CALLOC(bundle_count,seq->nsource_seq,int);
  LOOP (i,path_length) /* COUNT #POSITIONS OF EACH SEQ ARE IN path */
    for (source= &seq->letter[path[i]].source;source;source=source->more)
      bundle_count[source->iseq]++;

  LOOP (i,seq->nsource_seq) {/* FOR EACH SEQ OVER THRESHOLD, adjust weight*/
    /* printf("reweigh:\t%s\t%d/%d %d\n",seq->source_seq[i].name, */
      /* bundle_count[i],seq->source_seq[i].length,seq->source_seq[i].weight); */
    if (seq->source_seq[i].bundle_id<0 /* NOT YET BUNDLED */
	&& seq->source_seq[i].length*minimum_fraction <= bundle_count[i]) {
      /* printf("   +++++++++++++++++"); */
      seq->source_seq[i].weight = weight_multiplier*seq->source_seq[i].weight; /* Increase weight */
      nseq_in_bundle++;
      reweighted[i] = 1;
    } else {
	reweighted[i] = 0;
    }
    /* printf("\n"); */
  }

  FREE(bundle_count);
  return nseq_in_bundle; /* RETURN COUNT OF SEQUENCES IN BUNDLE */
}




/** assigns weights for bundling based upon /hb_weight arguments
    in source_seq titles */

void assign_hb_weights(int nsource_seq,LPOSourceInfo_T source_seq[])
{
  int i,weight;
  char *p;
  LOOP (i,nsource_seq) {
    if (source_seq[i].title &&
	(p=strstr(source_seq[i].title,"/hb_weight="))) {
      weight=atoi(p+11);
      if (weight!=0){ /* 0 COULD MEAN atoi FAILED TO PARSE ARG.  IGNORE IT*/
	source_seq[i].weight = weight;
	/* fprintf(stderr,"assigned weight=%d to %s\n",source_seq[i].weight,source_seq[i].name); */
      }
      else
	WARN_MSG(USERR,(ERRTXT,"hb_weight zero or unreadable: %s\nIgnored",p),"$Revision: 1.2 $");
    }
  }
}




/* Leena: multiplier for pulling out statistically dependent seqs */
//#define WEIGHT_MULTIPLIER 1000


/** generates the complete set of heaviest_bundle traversals of the the LPO
 seq, using iterative heaviest_bundle() and requiring that at least
 minimum_fraction of the positions in a sequence match the heaviest
 bundle path, for that sequence to be assigned to that bundle 
---------------------------------------------------------------
------------------------------------------------------------*/

/* ### NEW VERSION by Leena ###*/
void generate_lpo_bundles_MOD_LS(LPOSequence_T *seq,float minimum_fraction)
//void generate_lpo_bundles(LPOSequence_T *seq,float minimum_fraction)
{
  int nbundled=0,ibundle=0,path_length,iseq,count;
  LPOLetterRef_T *path=NULL;
  char name[256],title[1024];
  int *reweighted;  /* array for indicating which seqs have been reweighed */

  /* printf("#seq: %d\n", seq->nsource_seq); */
  /*  assign_hb_weights(seq->nsource_seq,seq->source_seq); TURN THIS ON!!*/
  while (nbundled < seq->nsource_seq) {/* PULL OUT BUNDLES ONE BY ONE */
    /* printf("more... %d %d\n", nbundled, seq->nsource_seq); */
    //FREE (path);
    //FREE (reweighted);

    /*Leena: Find first unbundled sequence and increase its weight*/
    int primary_seq = -1;
    for(int i = 0; i< seq->nsource_seq; i++) {
	if (seq->source_seq[i].bundle_id < 0) {
	    seq->source_seq[i].weight = WEIGHT_MULTIPLIER*seq->source_seq[i].weight;
	    primary_seq = i;
	    break;
	}
    }
    if (primary_seq < 0)
      goto premature_warning;
    //if (path)
    //  FREE (path);
    //  path=NULL;
    
    path=heaviest_bundle(seq->length,seq->letter,/*GET NEXT HEAVIEST BUNDLE*/
			 seq->nsource_seq,seq->source_seq,&path_length);
    /* ??!? FAILED TO FIND A BUNDLE ??? */
    if (!path || path_length<4){
      FREE (path);
      goto premature_warning;
    }
    /*Leena: reweigh the similar seqs*/
    //if (reweighted)  
    //  FREE (reweighted);
    //  reweighted=NULL;
    CALLOC(reweighted, seq->nsource_seq, int);
    /* Reset the weight of the primary seq so that it won't be multiplied twice */
    seq->source_seq[primary_seq].weight = seq->source_seq[primary_seq].weight/WEIGHT_MULTIPLIER;                                                //   THIS LINE GOT COMMENTED OUT FOR TESTING // AND NOW RE ENABLED
    /* Adjust weights */
    count=reweigh_sequences_in_bundle(path_length,path,seq,WEIGHT_MULTIPLIER, minimum_fraction, reweighted);
    //if (path)
    //  FREE (path);
    //  path=NULL;
    FREE (path);
    path=heaviest_bundle(seq->length,seq->letter,/*GET NEXT HEAVIEST BUNDLE*/
			 seq->nsource_seq,seq->source_seq,&path_length);
    /* ??!? FAILED TO FIND A BUNDLE ??? */
    if (!path || path_length<4){
      FREE (path);
      FREE (reweighted);
      goto premature_warning;
    }

    sprintf(name,"CONSENS%d",ibundle);
    /* NEXT, MARK SEQUENCES THAT FIT THIS BUNDLE ADEQUATELY */
    count=assign_sequence_bundle_id(path_length,path,seq,ibundle,
				    minimum_fraction);
    /*Leena: readjust the weights*/
    for(int i = 0; i< seq->nsource_seq; i++) {
	    if(reweighted[i])
	    seq->source_seq[i].weight = seq->source_seq[i].weight/WEIGHT_MULTIPLIER;
    }
    
    sprintf(title,"consensus produced by heaviest_bundle, containing %d seqs",count); /* DON'T INCLUDE CONSENSUS ITSELF IN THE COUNT! */
    iseq=add_path_sequence(path_length,path,seq,name,title);/*BUILD CONSENSUS*/
    seq->source_seq[iseq].bundle_id=ibundle++; /* INCREMENT BUNDLE ID */
    nbundled+=count+1; /* KEEP TRACK OF TOTAL SEQUENCES BUNDLED */

    FREE (reweighted);
    FREE (path);

    if (count<1) {
    premature_warning:
      //fprintf(stderr,"*** WARNING: bundling ended prematurely after %d bundles.\nNo sequences fit inside this last bundle.\nA total of %d sequences incuding consensus were bundled.\n\n",ibundle,nbundled);
      break;
    }
  }
  //if (reweighted)
  //  FREE (reweighted);
  //  reweighted=NULL;
  //if (path)
  //  FREE (path);
  //  path=NULL;
  
}

/* ### ORIGINAL VERSION ### */

void generate_lpo_bundles_ORIGINAL(LPOSequence_T *seq,float minimum_fraction)
{
  int nbundled=0,ibundle=0,path_length,iseq,count;
  LPOLetterRef_T *path=NULL;
  char name[256],title[1024];

  //printf("#seq: %d\n", seq->nsource_seq);

  /*  assign_hb_weights(seq->nsource_seq,seq->source_seq); TURN THIS ON!!*/
  while (nbundled < seq->nsource_seq) {/* PULL OUT BUNDLES ONE BY ONE */
    //printf("more... %d %d\n", nbundled, seq->nsource_seq);
    FREE (path);
    path=heaviest_bundle(seq->length,seq->letter,seq->nsource_seq,seq->source_seq,&path_length);/*GET NEXT HEAVIEST BUNDLE*/
    if (!path || path_length<10) /* ??!? FAILED TO FIND A BUNDLE ??? */
      goto premature_warning;
    sprintf(name,"CONSENS%d",ibundle);
    /* NEXT, MARK SEQUENCES THAT FIT THIS BUNDLE ADEQUATELY */
    count=assign_sequence_bundle_id(path_length,path,seq,ibundle,minimum_fraction);
    sprintf(title,"consensus produced by heaviest_bundle, containing %d seqs",count); /* DON'T INCLUDE CONSENSUS ITSELF IN THE COUNT! */
    iseq=add_path_sequence(path_length,path,seq,name,title);/*BUILD CONSENSUS*/
    seq->source_seq[iseq].bundle_id=ibundle++; /* INCREMENT BUNDLE ID */
    nbundled+=count+1; /* KEEP TRACK OF TOTAL SEQUENCES BUNDLED */

    if (count<1) {
    premature_warning:
      //fprintf(stderr,"*** WARNING: bundling ended prematurely after %d bundles.\nNo sequences fit inside this last bundle.\nA total of %d sequences incuding consensus were bundled.\n\n",ibundle,nbundled);
      break;
    }
  }
  FREE (path);
}



int compare_fun_small2big (const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}

int compare_fun_big2small (const void * a, const void * b) {
   return ( *(int*)b - *(int*)a );
}

/* Heaviest bundle algorithm with smart sequence weighting REWEIGH WHOLE BUNDLE*/
//void generate_lpo_bundles_SMART_REWEIGH_WHOLE_BUNDLE(LPOSequence_T *seq,float minimum_fraction)
void generate_lpo_bundles(LPOSequence_T *seq,float minimum_fraction)
{
  int nbundled=0,ibundle=0,path_length,iseq,count;
  LPOLetterRef_T *path=NULL;
  char name[256],title[1024];
  int *reweighted;  /* array for indicating which seqs have been reweighed */
  int graph_len = seq->length;
  int source_seqs = seq->nsource_seq;
  int *node_weights = NULL;
  int *node_weights_sorted = NULL;
  int *sequence_good_node_counts = NULL;
  CALLOC(node_weights, graph_len, int);
  CALLOC(node_weights_sorted, graph_len, int);
  CALLOC(sequence_good_node_counts, source_seqs, int);
  double top_node_percentage = 0.75;

  if (source_seqs <= 0){return;}

  int seq_count;
  int total_seq_len;
  int avg_seq_len;

  int absolute_minimum_good_node_score = 3;

  int ii;
  LOOP (ii,graph_len)
    node_weights[ii] = 0;

  LOOP (ii,source_seqs)
    sequence_good_node_counts[ii] = 0;

  /* printf("#seq: %d\n", seq->nsource_seq); */
  /*  assign_hb_weights(seq->nsource_seq,seq->source_seq); TURN THIS ON!!*/
  while (nbundled < seq->nsource_seq) {/* PULL OUT BUNDLES ONE BY ONE */
    /* printf("more... %d %d\n", nbundled, seq->nsource_seq); */
    //FREE (path);
    //FREE (reweighted);

    // Calculate average sequence length
    seq_count = 0;
    total_seq_len = 0;
    for (int gg = 0; gg < seq->nsource_seq; gg++){
      if (seq->source_seq[gg].weight > 0){
        seq_count += 1;
        total_seq_len += seq->source_seq[gg].length;
      }
    }

    if (seq_count <= 0){avg_seq_len = 0;}
    else{avg_seq_len = (int)round(total_seq_len/seq_count);}
    

    /*Leena: Find first unbundled sequence and increase its weight*/
    /*Miika: Find the "best" unbundled sequence and increase its weight*/
   
    // Reset node weights
    memset(node_weights, 0, graph_len*sizeof(int));
    memset(node_weights_sorted, 0, graph_len*sizeof(int));
    int weight = 0;
    int nonzeronodes = 0;
    int source_seq_id;
    // Recalculate node weights using updated sequence weights
    for (int a = 0; a < graph_len; a++){
      weight = 0;
      //printf("node %d has seqs ", a);
      LPOLetterSource_T *source_seq = &seq->letter[a].source;
      source_seq_id = source_seq->iseq;
      //printf("%d ", source_seq_id);
     
      if (seq->source_seq[source_seq_id].weight > 0){weight+=1;}

      while (source_seq=source_seq->more){
        source_seq_id = source_seq->iseq;
        //printf("%d ", source_seq_id);
        if (seq->source_seq[source_seq_id].weight > 0){weight+=1;}
      }

      node_weights[a] = weight;
      node_weights_sorted[a] = weight;
      if (weight > 0){nonzeronodes+=1;}
      //printf("and has weight %d \n", weight);
      
    }

    /* Figure out what the minimum good node score should be*/
    qsort(node_weights_sorted, graph_len, sizeof(int), compare_fun_big2small);
    if (avg_seq_len >= graph_len){avg_seq_len = graph_len-1;}
    int good_node_min_score = node_weights_sorted[avg_seq_len];

    if (good_node_min_score < absolute_minimum_good_node_score){good_node_min_score = absolute_minimum_good_node_score;}
    
    int cut = 0;
    int bonus = 0;
    for (int d = avg_seq_len-1; d >= 0; d--){
      if (node_weights_sorted[d] != good_node_min_score){break;}
      else {cut += 1;}
    }
    for (int e = avg_seq_len+1; e < seq->nsource_seq; e++){
      if (node_weights_sorted[e] != good_node_min_score){break;}
      else {bonus += 1;}
    }
    if (bonus >= cut){good_node_min_score = good_node_min_score + 1;}
    



    //printf("Lowest position 0: %f\n", top_node_percentage);
    //double bottom_percentage = 1 - top_node_percentage;
    //printf("Lowest position 1: %f\n", bottom_percentage);
    //double lowest_position_float = bottom_percentage * nonzeronodes;
    //printf("Lowest position 2: %f\n", lowest_position_float);
    //int lowest_position = (int)lowest_position_float + 1;
    //printf("Lowest position 3: %d\n", lowest_position);

    //int lowest_accepted_weight_position = (graph_len - nonzeronodes) + lowest_position;
    //if (lowest_accepted_weight_position < 0){lowest_accepted_weight_position = 0;}
    //if (lowest_accepted_weight_position >= graph_len){lowest_accepted_weight_position = graph_len-1;}
    //int lowest_accepted_weight = node_weights_sorted[lowest_accepted_weight_position];
    //printf("Lowest accepted weight is %d at lowest accepted position %d\n", lowest_accepted_weight, lowest_accepted_weight_position);

    //if (lowest_accepted_weight < 3){
    //  lowest_accepted_weight = 3;
      //printf("Lowest accepted weight raised to 3\n");
    //}


    int good_nodes_count = 0;
    memset(sequence_good_node_counts, 0, source_seqs*sizeof(int));

    // Calculate how many good nodes each sequence has
    for (int b = 0; b < graph_len; b++){
      //printf("%d -- %d\n", node_weights[b], lowest_accepted_weight);
      if (node_weights[b] < good_node_min_score){continue;}
      //printf("In good node %d following seqs: ", b);
      good_nodes_count += 1;
      LPOLetterSource_T *source_seq = &seq->letter[b].source;
      source_seq_id = source_seq->iseq;
      //printf("%d ", source_seq_id);
      if ((seq->source_seq[source_seq_id].weight > 0)&&(source_seq_id < source_seqs)){
        sequence_good_node_counts[source_seq_id]+=1;
        //printf("%d ", source_seq_id);
      }
      while (source_seq=source_seq->more){
        source_seq_id = source_seq->iseq;
        if ((seq->source_seq[source_seq_id].weight > 0)&&(source_seq_id < source_seqs)){
          //printf("%d ", source_seq_id);
          sequence_good_node_counts[source_seq_id]+=1;
        }
      }
      //printf("\n");
    }


    // Find the best sequence
    int best_seq_id = -1;
    int best_seq_score = 999999999;
    for (int c = 0; c < source_seqs; c++){
      if (sequence_good_node_counts[c] == 0){continue;}
      // Good nodes not in the sequence
      int current_score = good_nodes_count - sequence_good_node_counts[c];
      // Add nodes in the sequence that are not good
      current_score = current_score + (seq->source_seq[c].length - sequence_good_node_counts[c]);
      if(current_score < best_seq_score){
        best_seq_id = c;
        best_seq_score = current_score;
      }
      //printf("Sequence %d length is %d and has %d good nodes and score is: %d\n", c, seq->source_seq[c].length, sequence_good_node_counts[c], current_score);
    }

    //printf("Good node threshold is: %d\n", good_node_min_score);
    //printf("Good nodes count: %d\n", good_nodes_count);
    //printf("Best sequence is: %d\n", best_seq_id);

    int primary_seq = best_seq_id;

    if (primary_seq < 0){
      //FREE (node_weights);
      //FREE (sequence_good_node_counts);
      goto premature_warning;
    }

    seq->source_seq[primary_seq].weight = WEIGHT_MULTIPLIER*seq->source_seq[primary_seq].weight;                             // DIS NEW STUFF

    /*
    int primary_seq = -1;
    for(int i = 0; i< seq->nsource_seq; i++) {
      if (seq->source_seq[i].bundle_id < 0) {
        seq->source_seq[i].weight = WEIGHT_MULTIPLIER*seq->source_seq[i].weight;
        primary_seq = i;
        break;
      }
    }
    */
    
    //if (path)
    //  FREE (path);
    //  path=NULL;
    
    path=heaviest_bundle(seq->length,seq->letter,/*GET NEXT HEAVIEST BUNDLE*/
       seq->nsource_seq,seq->source_seq,&path_length);
    /* ??!? FAILED TO FIND A BUNDLE ??? */
    if (!path || path_length<10){
      FREE (path);
      //FREE (node_weights);
      //FREE (sequence_good_node_counts);
      goto premature_warning;
    }
    /*Leena: reweigh the similar seqs*/
    //if (reweighted)  
    //  FREE (reweighted);
    //  reweighted=NULL;
    CALLOC(reweighted, seq->nsource_seq, int);
    /* Reset the weight of the primary seq so that it won't be multiplied twice */


    /*
    printf("BEFORE REWEIGHING\n");
    for (int pop = 0; pop < source_seqs; pop++){
      printf("SEQ %d WEIGHT = %d\n", pop, seq->source_seq[pop].weight);
    }*/

    seq->source_seq[primary_seq].weight = seq->source_seq[primary_seq].weight / WEIGHT_MULTIPLIER;                                               // ENABLE THIS SINCE IT LOOKS BAD WITHOUT IT

    /*
    printf("AFTER BEST SEQ WEIGHT LOWERED\n");
    for (int pop = 0; pop < source_seqs; pop++){
      printf("SEQ %d WEIGHT = %d\n", pop, seq->source_seq[pop].weight);
    }*/

    /* Adjust weights */
    count=reweigh_sequences_in_bundle(path_length,path,seq,WEIGHT_MULTIPLIER, minimum_fraction, reweighted);

    /*
    printf("AFTER BUNDLE WEIGHT INCREASED\n");
    for (int pop = 0; pop < source_seqs; pop++){
      printf("SEQ %d WEIGHT = %d\n", pop, seq->source_seq[pop].weight);
    }*/

    //if (path)
    //  FREE (path);
    //  path=NULL;
    FREE (path);
    path=heaviest_bundle(seq->length,seq->letter,/*GET NEXT HEAVIEST BUNDLE*/
       seq->nsource_seq,seq->source_seq,&path_length);
    /* ??!? FAILED TO FIND A BUNDLE ??? */
    if (!path || path_length<10){
      FREE (path);
      FREE (reweighted);
      //FREE (node_weights);
      //FREE (sequence_good_node_counts);
      goto premature_warning;
    }

    sprintf(name,"CONSENS%d",ibundle);
    /* NEXT, MARK SEQUENCES THAT FIT THIS BUNDLE ADEQUATELY */
    count=assign_sequence_bundle_id(path_length,path,seq,ibundle,
            minimum_fraction);
    /*Leena: readjust the weights*/
    for(int i = 0; i< seq->nsource_seq; i++) {
      if(reweighted[i])
      seq->source_seq[i].weight = seq->source_seq[i].weight/WEIGHT_MULTIPLIER;
    }
    
    sprintf(title,"consensus produced by heaviest_bundle, containing %d seqs",count); /* DON'T INCLUDE CONSENSUS ITSELF IN THE COUNT! */
    iseq=add_path_sequence(path_length,path,seq,name,title);/*BUILD CONSENSUS*/
    seq->source_seq[iseq].bundle_id=ibundle++; /* INCREMENT BUNDLE ID */
    nbundled+=count+1; /* KEEP TRACK OF TOTAL SEQUENCES BUNDLED */

    FREE (reweighted);
    FREE (path);
    //FREE (node_weights);
    //FREE (sequence_good_node_counts);

    if (count<1) {
    premature_warning:
      //fprintf(stderr,"*** WARNING: bundling ended prematurely after %d bundles.\nNo sequences fit inside this last bundle.\nA total of %d sequences incuding consensus were bundled.\n\n",ibundle,nbundled);
      break;
    }
  }
  FREE (node_weights);
  FREE (node_weights_sorted);
  FREE (sequence_good_node_counts);
}