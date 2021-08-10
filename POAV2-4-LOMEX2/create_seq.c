

#include "default.h"
#include "seq_util.h"




void save_sequence_fields(Sequence_T *seq,
			  char seq_name[],char seq_title[],int length)
{
  STRNCPY(seq->name,seq_name,SEQUENCE_NAME_MAX);
  if (seq_title)
    seq->title=strdup(seq_title);
  else
    seq->title=strdup("untitled");
  seq->length=length; /* SAVE LENGTH */
}



int create_seq(int nseq,Sequence_T **p_seq,
	       char seq_name[],char seq_title[],char tmp_seq[],
	       int do_switch_case)
{
  //printf("begin seq make\n");
  int i,j;
  Sequence_T *seq;

  //printf("sm1\n");

  REBUFF(*p_seq,nseq,SEQUENCE_BUFFER_CHUNK,Sequence_T); /* ALLOCATE MEMORY*/
  //printf("sm2\n");
  seq= (*p_seq)+nseq; /* SET POINTER TO NEWLY ALLOCATED ELEMENT */
  //printf("sm3\n");
  
  for (i=j=0;tmp_seq[i];i++) /* ELIMINATE WHITE SPACE */
    if (!isspace(tmp_seq[i]))
      tmp_seq[j++]=tmp_seq[i];
  //printf("sm4\n");  
  tmp_seq[j]='\0'; /* TERMINATE COMPRESSED STRING*/
  //printf("sm5\n");
  seq->sequence=strdup(tmp_seq); /* SAVE A DYNAMIC COPY */
  //printf("sm6\n");
  save_sequence_fields(seq,seq_name,seq_title,j);
  //printf("sm7\n");
  switch (do_switch_case) {
  case switch_case_to_lower:
    LOOP (i,seq->length)
      seq->sequence[i]=tolower(tmp_seq[i]);
    break;
  case switch_case_to_upper:
    LOOP (i,seq->length)
      seq->sequence[i]=toupper(tmp_seq[i]);
    break;
  }

  //printf("end seq make\n");
  return 1;
}


