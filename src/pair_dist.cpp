#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <math.h>
#include "pair_dist.h"

using namespace std;

#define FOUR_THIRDS (4.0/3.0)
Pair_dist::Pair_dist()
{
  cerr<<"Called to Pair_dist default constructor\n";
}


Pair_dist::Pair_dist(Sequence_dataset *the_seq, Exchange *cexchange)
{
  int i, j;
  curr_sequence=the_seq;
  curr_exchange=cexchange;

	ambigs_ignored=FALSE;
  total_remaining=curr_exchange->get_num_taxa();

  pos=-1;
  joined=new int [curr_exchange->get_num_taxa()];


  dists= new distances [curr_exchange->get_num_taxa()];  
  for(i=0; i<curr_exchange->get_num_taxa(); i++)
    {
      dists[i].id=i;
      dists[i].num_com=1;
      dists[i].t_distance= new double [curr_exchange->get_num_taxa()]; 
      for (j=0; j<curr_exchange->get_num_taxa(); j++)	  
	dists[i].t_distance[j]=0;  
    }

} //End Pair_dist::Pair_dist



void Pair_dist::reset_seqs(Sequence_dataset *new_seq)
{
	if (new_seq->Num_sequences() != curr_exchange->get_num_taxa()) cerr<<"Cannot change number of taxa when switching alignments in Pair_dist\n";
	else curr_sequence=new_seq;
}

void Pair_dist::find_dist ()
{
  int i, m;

  for(i=0; i<curr_sequence->Num_sequences(); i++)      
    for(m=i+1; m<curr_sequence->Num_sequences(); m++)
	pairwise_dist(i,m);
	
}//end Pair_dist::find_dist




double Pair_dist::return_dist(int taxa_1, int taxa_2)
{
  int i, j;
  double retval;
  
  for(i=0; i<total_remaining; i++)
    {
      if (dists[i].id==taxa_1)
	{
	  for (j=0; j<total_remaining; j++)
	    if (dists[j].id==taxa_2)
	      retval=dists[i].t_distance[j];
	}
      
    }
  return(retval);
} //End Pair_dist::return_dist




double Pair_dist::return_closest_pair(int &taxa_1, int &taxa_2, BOOL &joined_1, BOOL &joined_2)
{
  int i, j;
  double best;

  best=dists[0].t_distance[1];
  taxa_1=dists[0].id;
  taxa_2=dists[1].id;

  for(i=0; i<total_remaining; i++)
    for(j=2; j<total_remaining; j++)
      {
	if(i!=j)
	  if (dists[i].t_distance[j]<best)
	    {
	      best=dists[i].t_distance[j];
	      taxa_1=dists[i].id;
	      taxa_2=dists[j].id;
	    }

      }

  joined_1=is_joined(taxa_1);
  joined_2=is_joined(taxa_2);

  return(best);
} //End Pair_dist::return_closest_pair




void Pair_dist::collapse_matrix (int on_1, int on_2)
{
  int i, j, temp, loc_1, loc_2;
 
  if (on_2<on_1)
    {
      on_2=temp;
      on_2=on_1;
      on_1=temp;
    }
 
  i=0;
  while(dists[i].id!=on_1)
    i++;
  
  loc_1=i;

  
  while(dists[i].id!=on_2)
    i++;
  
  loc_2=i;

  for (i=0; i<total_remaining; i++)
      if (dists[i].id!=on_1)
	  dists[loc_1].t_distance[i]=dists[i].t_distance[loc_1]=
	    ((dists[loc_1].num_com*dists[loc_1].t_distance[i])+
	     (dists[loc_2].num_com*dists[loc_2].t_distance[i]))/
	    (dists[loc_1].num_com+dists[loc_2].num_com);

  for (i=0; i<loc_2; i++)      
      for(j=loc_2; j<total_remaining; j++)
	dists[i].t_distance[j]=dists[i].t_distance[j+1];

  for (i=loc_2; i<total_remaining-1; i++)
    {
      dists[i].id=dists[i+1].id;
      for(j=0; j<loc_2; j++)
	dists[i].t_distance[j]=dists[i+1].t_distance[j];
      for(j=loc_2; j<total_remaining; j++)
	dists[i].t_distance[j]=dists[i+1].t_distance[j+1];
      
    }
  pos++;
  joined[pos]=on_1;
  pos++;
  joined[pos]=on_2;

  total_remaining--;
}




void Pair_dist::print_dists()
{
  int i, j, m;
  char name_copy[700];

  for(i=0; i<total_remaining; i++)
       {
	 for(m=0; m<10; m++)
	   {
	     strcpy(name_copy, (*curr_sequence)[dists[i].id].Sequence_name());
	     if (strlen(name_copy)>m)
	       cout<<name_copy[m];
	     else
	       cout<<" ";
	   }
	 cout<<"  ";
	 
	 for(m=0; m<total_remaining; m++)
	   {
	     if (i==m)
	       cout<<"    --    ";
	     else
	       cout<<setw(8)<<setprecision(4)
		   <<(1.0*dists[i].t_distance[m])<<"  ";
	   }
	 cout<<endl;
	 
       }
  cout<<endl;
  
} //End Pair_dist::print_dists



Pair_dist::~Pair_dist()
{
  int i;
  
	//cout<<"Pair_dist destructor\n"<<flush;

	
  for(i=0; i<curr_sequence->Num_sequences(); i++)
    delete[] dists[i].t_distance; 
  
  delete[] dists;
  delete[] joined;
}



BOOL Pair_dist::is_joined(int taxa)
{
  int i;
  BOOL retval=FALSE;

  for(i=0; i<=pos; i++)
    if(joined[i]==taxa)
      retval=TRUE;
  return(retval);

}




//Nucleotide_dist functions
void Nucleotide_dist::allocate_diff_matrix()
{
  int i,j;  
  num_diffs=new int* [curr_exchange->get_num_taxa()];
  
  for(i=0; i<curr_exchange->get_num_taxa(); i++)
    {
      num_diffs[i]=new int [curr_exchange->get_num_taxa()];
  
      for (j=0; j<curr_exchange->get_num_taxa(); j++) 
	num_diffs[i][j]=0;
	
       }
} //End Nucleotide_dist::allocate_diff_matrix




Nucleotide_dist::~Nucleotide_dist()
{
  int i;
  
	//cout<<"Nucleotide_dist destructor\n"<<flush;
	
  for (i=0; i<curr_exchange->get_num_taxa(); i++)
    delete[] num_diffs[i];
  delete[] num_diffs;
}




void Nucleotide_dist::handle_gaps(int taxa, int character, int &start, int &end)
{
  if (num_to_base((*curr_sequence)[taxa][character])=='-')
    {
      start=0; 
      end=4;
    }
  else
    start=end=(*curr_sequence)[taxa][character];
} //End Nucleotide_dist::handle_gaps




//Codon_dist functions
Codon_dist::Codon_dist(Sequence_dataset *the_seq, Exchange *cexchange) : Pair_dist(the_seq, cexchange)
{
  allocate_diff_matrix();
  curr_gen_code=create_genetic_code(cexchange);
}


Codon_dist::~Codon_dist()
{
  int i;
  
  for(i=0; i<curr_sequence->Num_sequences(); i++)
      delete[] total_diffs[i];

  delete[] total_diffs;
  delete curr_gen_code;
}



void Codon_dist::allocate_diff_matrix()
{
  int i,j;  
  total_diffs=new double* [curr_exchange->get_num_taxa()];
  
  for(i=0; i<curr_exchange->get_num_taxa(); i++)
    {
      total_diffs[i]=new double [curr_exchange->get_num_taxa()];
  
      for (j=0; j<curr_exchange->get_num_taxa(); j++) 
	     total_diffs[i][j]=0.0;
       }

} //End Codon_dist::allocate_diff_matrix




int Codon_dist::num_inters_wo_stops(int codon1[3], int codon2[3], int num_diffs)
{
  int i, retval=0, inter_codon1[3], inter_codon2[3], dum_diff1, dum_diff2, dum_diff3;

  switch (num_diffs)
    {
    case 2:
      for (i=0; i<2; i++)
	{
	  intermediate_codon(codon1, codon2, num_diffs, i, dum_diff1, dum_diff2, dum_diff3,
			  inter_codon1, inter_codon2);
      
	  if(curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon1))==FALSE)
	    retval++;
	}
      break;
    case 3:
      for (i=0; i<6; i++)
	{
	  intermediate_codon(codon1, codon2, num_diffs, i, dum_diff1, dum_diff2, dum_diff3,
			  inter_codon1, inter_codon2);
      
	  if(curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon1))==FALSE
	     && curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon2))==FALSE)
	    retval++;
	}
      break;
    }
  
  return(2*retval);
}  //End Codon_dist::num_inters_wo_stops




void Codon_dist::intermediate_codon(int codon1[3], int codon2[3], int num_diffs, int intermed_num, int &diff1,
				    int &diff2, int &diff3, int intermed_codon[3], int intermed_codon1[3])
{
  int i, temp;
  BOOL diffs[3];

  curr_gen_code->diff_pos(codon1, codon2, diffs);

  switch (num_diffs)
    {
    case 2:
      diff1=0;
      while(diffs[diff1]==FALSE)
	diff1++;
      diff2=diff1+1;
      while(diffs[diff2]==FALSE)
	diff2++;
      if (intermed_num!=1)
	{
	  temp=diff1;
	  diff1=diff2;
	  diff2=temp;
	}

      for (i=0; i<diff1; i++)
	intermed_codon[i]=codon1[i];

      intermed_codon[diff1]=codon2[diff1];

      for (i=diff1+1; i<3; i++)
	intermed_codon[i]=codon1[i];
      break;
    case 3:
      switch (intermed_num)
	{
	case 0:
	  diff1=0;
	  diff2=1;
	  diff3=2;
	  break;
	case 1:
	  diff1=0;
	  diff2=2;
	  diff3=1;
	  break;
	case 2:
	  diff1=1;
	  diff2=0;
	  diff3=2;
	  break;
	case 3:
	  diff1=1;
	  diff2=2;
	  diff3=0;
	  break;
	case 4:
	  diff1=2;
	  diff2=0;
	  diff3=1;
	  break;
	case 5:
	  diff1=2;
	  diff2=1;
	  diff3=0;
	  break;
	}
      
      for (i=0; i<diff1; i++)
	intermed_codon[i]=codon1[i];

      intermed_codon[diff1]=codon2[diff1];

      for (i=diff1+1; i<3; i++)
	intermed_codon[i]=codon1[i];

      for (i=0; i<diff2; i++)
	intermed_codon1[i]=intermed_codon[i];

      intermed_codon1[diff2]=codon2[diff2];

      for (i=diff2+1; i<3; i++)
	intermed_codon1[i]=intermed_codon[i];
      break;
    }

}  //End Codon_dist::intermediate_codon




void Codon_dist::gap_possible_codons(int taxa, int character, int &numgaps, int &num_gap_codons, int &gap1, int &gap2)
{
  int i, j, k, gaps=1, codon[3];
  BOOL gap_pos[3];
  
  numgaps=0;
  num_gap_codons=1;

  for (i=0; i<3; i++)
    {
      codon[i]=(*curr_sequence)[taxa][3*character+i];
	if (num_to_base((*curr_sequence)[taxa][3*character+i])=='-')
	  {
	    gap_pos[i]=TRUE;
	    numgaps++;
	  }
	else
	  gap_pos[i]=FALSE;
    }

  if(numgaps>0)
    {
      num_gap_codons=4;
      i=0;
      while(gap_pos[i]==FALSE)
	i++;

      gap1=i;
      if(numgaps>1)
	{
	  num_gap_codons=16;
	  i++;
	  while(gap_pos[i]==FALSE)
	    i++;
	  gap2=i;
	  if(numgaps==3)
	    num_gap_codons=64;
	}
    }

  //Stop codon corrections:
  switch (numgaps)
    {
    case 1:
      for (i=0; i<4; i++)
	{
	  codon[gap1]=i;
	  if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon))==TRUE)
	    num_gap_codons--;

	}
      break;
    case 2:
      for (i=0; i<4; i++)
	{
	  codon[gap1]=i;
	  for (j=0; j<4; j++)
	    {
	      codon[gap2]=j;
	      if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon))==TRUE)
		num_gap_codons--;
	    }
	}
      break;
     
    case 3:
      for (i=0; i<4; i++)
	{
	  codon[0]=i;
	  for (j=0; j<4; j++)
	    {
	      codon[1]=j;
	      for (k=0; k<4; k++)
		{
		  codon[2]=k;
		  if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon))==TRUE)
		    num_gap_codons--;
		}
	    }
	}
      break;
      
    }
} //End gap_possible_codons




void Codon_dist::get_gap_codon(int codon[3], int gcodonnum, int numgaps, int gap1, int gap2)
{
  switch(numgaps)
    {
    case 1:
      codon[gap1]=gcodonnum;
      break;
    case 2:
      codon[gap1]=gcodonnum/4;
      codon[gap2]=gcodonnum%4;
      break;
    case 3:
      codon[0]=gcodonnum/16;
      codon[1]=(gcodonnum-16*(gcodonnum/16))/4;
      codon[2]=(gcodonnum-16*(gcodonnum/16))%4;
      break;
    }
} //End Codon_dist::get_gap_codon





void Amino_Acid_dist::handle_gaps(int taxa, int character, int &start, int &end)
{
  if (num_to_aa((*curr_sequence)[taxa][character])=='-')
    {
      start=0; 
      end=20;
    }
  else
    start=end=(*curr_sequence)[taxa][character];
} //End Amino_Acid_dist::handle_gaps
  



void Amino_Acid_dist::allocate_diff_matrix()
{
  int i,j;  
  num_diffs=new int* [curr_exchange->get_num_taxa()];
  
  for(i=0; i<curr_exchange->get_num_taxa(); i++)
    {
      num_diffs[i]=new int [curr_exchange->get_num_taxa()];
      
      for (j=0; j<curr_exchange->get_num_taxa(); j++) 
	num_diffs[i][j]=0;
      
    }
}

Percent_dist::~Percent_dist()
{
	//cout<<"Percent_dist destructor\n"<<flush;	
}

void Percent_dist::pairwise_dist(int taxa_1, int taxa_2)
{
 int i, num_ignore=0, size;
   for(i=0; i<(*curr_sequence)[0].Sequence_size(); i++)
     {
		 if ((ambigs_ignored == FALSE) || 
			 ((base_is_ambig((*curr_sequence)[taxa_1][i]) ==FALSE) && (base_is_ambig((*curr_sequence)[taxa_2][i]) ==FALSE)) ) {	 
			if((*curr_sequence)[taxa_1][i] != 
				(*curr_sequence)[taxa_2][i])
				{
					num_diffs[taxa_1][taxa_2]++;
					num_diffs[taxa_2][taxa_1]++;
				}
		 }
       if ((base_is_ambig((*curr_sequence)[taxa_1][i]) ==TRUE) || 
		   (base_is_ambig((*curr_sequence)[taxa_2][i]) ==TRUE)) num_ignore++;
     }

	size= (*curr_sequence)[0].Sequence_size();
	
	if (ambigs_ignored == TRUE)
		size -= num_ignore;
	
   dists[taxa_1].t_distance[taxa_2]=(1.0*num_diffs[taxa_1][taxa_2])/size;
	dists[taxa_2].t_distance[taxa_1]=(1.0*num_diffs[taxa_1][taxa_2])/size;    
}


JC_dist::JC_dist() : Nucleotide_dist() {}

JC_dist::JC_dist(Sequence_dataset *the_seq, Exchange *cexchange) : Nucleotide_dist(the_seq, cexchange) {}



void JC_dist::pairwise_dist(int taxa_1, int taxa_2)
{
 int i;

 for(i=0; i<(*curr_sequence)[0].Sequence_size(); i++)
   {
     if((*curr_sequence)[taxa_1][i] != 
	(*curr_sequence)[taxa_2][i])
       {
	 num_diffs[taxa_1][taxa_2]++;
	 num_diffs[taxa_2][taxa_1]++;
       }
     
   }
 
 dists[taxa_1].t_distance[taxa_2]=dists[taxa_2].t_distance[taxa_1]=
   -0.75*log(1-(FOUR_THIRDS*((1.0*num_diffs[taxa_1][taxa_2])/
    (*curr_sequence)[0].Sequence_size())));

}


NG_86_Synon_dist::NG_86_Synon_dist(): Codon_dist()
{
}

NG_86_Synon_dist::NG_86_Synon_dist(Sequence_dataset *the_seq, Exchange *cexchange) : Codon_dist(the_seq, cexchange)
{
} 




void NG_86_Synon_dist::pairwise_dist(int taxa_1, int taxa_2)
{
 int i, j, k, codon1[3], codon2[3], numgap1, numgap2, gapcodons1, gapcodons2, gap1a, gap2a, gap1b, gap2b;
 double diffs;
 BOOL gaps1[3], gaps2[3];

 possible_diffs=0.0;
 
 total_diffs[taxa_1][taxa_2]=total_diffs[taxa_2][taxa_1]=0.0;

 for(i=0; i<((*curr_sequence)[0].Sequence_size()/3); i++)
   {
     gap_possible_codons(taxa_1, i, numgap1, gapcodons1, gap1a, gap2a);
     gap_possible_codons(taxa_2, i, numgap2, gapcodons2, gap1b, gap2b);

     for(j=0; j<3; j++)
       {
	 codon1[j]=(*curr_sequence)[taxa_1][3*i+j];
	 codon2[j]=(*curr_sequence)[taxa_2][3*i+j];
       }
     for (j=0; j<gapcodons1; j++)
       {
	 get_gap_codon(codon1, j, numgap1, gap1a, gap2a);
	 if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon1))==FALSE)
	   {
	     for (k=0; k<gapcodons2; k++)
	       {
		 get_gap_codon(codon2, k, numgap2, gap1b, gap2b);
		 
		 if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon2))==FALSE)
		  { 
		    possible_diffs+=(1.0*curr_gen_code->possible_synon_subs(codon1)/
						       (gapcodons1*gapcodons2));
		    possible_diffs+=(1.0*curr_gen_code->possible_synon_subs(codon2)/
						       (gapcodons1*gapcodons2));
		    
		    diffs=(1.0*curr_gen_code->synon_paths(curr_gen_code->multiple_subs(codon1, codon2),
									    codon1, codon2)/(gapcodons1*gapcodons2));
		    
		    total_diffs[taxa_1][taxa_2]+=diffs;
		    total_diffs[taxa_2][taxa_1]+=diffs;
		  }
	       }
	   }
       }
   }

 possible_diffs=possible_diffs/2.0;

// cout<<"Nei Synonymous dist: total_diffs: "<<total_diffs[taxa_1][taxa_2]<<" Possible: "<<possible_diffs<<endl;
 
 if (total_diffs[taxa_1][taxa_2]!=0)
   dists[taxa_1].t_distance[taxa_2]=dists[taxa_2].t_distance[taxa_1]=
     -0.75*log(1.0-FOUR_THIRDS*(total_diffs[taxa_1][taxa_2]/
			       possible_diffs));
 else
   dists[taxa_1].t_distance[taxa_2]=dists[taxa_2].t_distance[taxa_1]=0;
 
} //End NG_86_Synon_dist::pairwise_dist



NG_86_Non_Synon_dist::NG_86_Non_Synon_dist(): Codon_dist()
{
}

NG_86_Non_Synon_dist::NG_86_Non_Synon_dist(Sequence_dataset *the_seq, Exchange *cexchange) : 
Codon_dist(the_seq, cexchange)
{
} 



void NG_86_Non_Synon_dist::pairwise_dist(int taxa_1, int taxa_2)
{
 int i, j, k, codon1[3], codon2[3], numgap1, numgap2, gapcodons1, gapcodons2, gap1a, gap2a, gap1b, gap2b;
 double diffs;

 possible_diffs=0.0;

 total_diffs[taxa_1][taxa_2]=total_diffs[taxa_2][taxa_1]=0.0;
 
 for(i=0; i<((*curr_sequence)[0].Sequence_size()/3); i++)
   {
       gap_possible_codons(taxa_1, i, numgap1, gapcodons1, gap1a, gap2a);
       gap_possible_codons(taxa_2, i, numgap2, gapcodons2, gap1b, gap2b);
    
     for(j=0; j<3; j++)
       {
	 codon1[j]=(*curr_sequence)[taxa_1][3*i+j];
	 codon2[j]=(*curr_sequence)[taxa_2][3*i+j];
       }
    
    for (j=0; j<gapcodons1; j++)
       {
	 get_gap_codon(codon1, j, numgap1, gap1a, gap2a);

	 if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon1))==FALSE)
	   {
	     for (k=0; k<gapcodons2; k++)
	       {
		 get_gap_codon(codon2, k, numgap2, gap1b, gap2b);
		 
		 if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon2))==FALSE)
		   {
		     possible_diffs+=(1.0*curr_gen_code->possible_non_synon_subs(codon1)/
							    (gapcodons1*gapcodons2));
		     possible_diffs+=(1.0*curr_gen_code->possible_non_synon_subs(codon2)/
							    (gapcodons1*gapcodons2));
		     
		     diffs=(1.0*curr_gen_code->non_synon_paths(
							       curr_gen_code->multiple_subs
							       (codon1, codon2), codon1, codon2)/
			    (gapcodons1*gapcodons2));
		     
		     total_diffs[taxa_1][taxa_2]+=diffs;
		     total_diffs[taxa_2][taxa_1]+=diffs;
		   }
	       }
	   }
       }
   }
 
 possible_diffs=possible_diffs/2.0;


 if (total_diffs[taxa_1][taxa_2]!=0)
   dists[taxa_1].t_distance[taxa_2]=dists[taxa_2].t_distance[taxa_1]=
     -0.75*log(1-(FOUR_THIRDS*((total_diffs[taxa_1][taxa_2])/
			       possible_diffs)));
 else
   dists[taxa_1].t_distance[taxa_2]=dists[taxa_2].t_distance[taxa_1]=0;
 
} //End NG_86_Non_Synon_dist::pairwise_dist




L_93_Synon_dist::L_93_Synon_dist(): Codon_dist()
{
}

L_93_Synon_dist::L_93_Synon_dist(Sequence_dataset *the_seq, Exchange *cexchange): Codon_dist(the_seq, cexchange)
{
}


void L_93_Synon_dist::pairwise_dist(int taxa_1, int taxa_2)
{
  int i, j, k,l, codon1[3], codon2[3], inter_codon[3], inter_codon2[3], num_diffs, diff1, diff2, 
    diff3, curr_degen,  numgap1, numgap2, gapcodons1, gapcodons2, gap1a, gap2a, gap1b, gap2b; 
  double l_level[5], s_level[5], v_level[5], temp_s, temp_v, valid_steps;


  for(i=0; i<5; i++)
    l_level[i]=s_level[i]=v_level[i]=0;


  for(i=0; i<((*curr_sequence)[0].Sequence_size()/3); i++)
    {
      gap_possible_codons(taxa_1, i, numgap1, gapcodons1, gap1a, gap2a);
      gap_possible_codons(taxa_2, i, numgap2, gapcodons2, gap1b, gap2b);
    

      for(j=0; j<3; j++)
	{
	  codon1[j]=(*curr_sequence)[taxa_1][3*i+j];
	  codon2[j]=(*curr_sequence)[taxa_2][3*i+j];
	}

      for (j=0; j<gapcodons1; j++)
	{
	  get_gap_codon(codon1, j, numgap1, gap1a, gap2a);

	  if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon1))==FALSE)
	   {      
	     for (k=0; k<gapcodons2; k++)
	       {
		 get_gap_codon(codon2, k, numgap2, gap1b, gap2b);
		 
		 if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon2))==FALSE)
		   {
		     num_diffs=curr_gen_code->multiple_subs(codon1, codon2);
		     
		     for(l=0; l<3; l++)
		       {
			 l_level[curr_gen_code->Li_93_degen(codon1, l)]+=1.0/(gapcodons1*gapcodons2);
			 l_level[curr_gen_code->Li_93_degen(codon2, l)]+=1.0/(gapcodons1*gapcodons2);
		       }
		     
		     switch (num_diffs)
		       {
		       case 1:
			 //cout<<"Codon 1: "<<codon1[0]<<codon1[1]<<codon1[2]<<" Codon 2: "<<codon2[0]<<codon2[1]<<codon2[2]
			 //  <<" Diffs: "<<num_diffs;
			 diff1=0;
			 while(codon1[diff1]==codon2[diff1])
			   diff1++;
			 
			 curr_degen=curr_gen_code->Li_93_degen(codon1, diff1);			 
			 temp_s=temp_v=0;

			 curr_gen_code->count_by_degen(codon1, codon2, diff1, temp_s, temp_v);			 
			 s_level[curr_degen]+=temp_s/(2.0*gapcodons1*gapcodons2);
			 v_level[curr_degen]+=temp_v/(2.0*gapcodons1*gapcodons2);
			 
			 // cout<<" Degen: "<<curr_degen<<" Forward S cnt: "<<temp_s<<" Forward V cnt: "<<temp_v;
			 
			 curr_degen=curr_gen_code->Li_93_degen(codon2, diff1);
			 temp_s=temp_v=0;
			 
			 curr_gen_code->count_by_degen(codon2, codon1, diff1, temp_s, temp_v);
			 s_level[curr_degen]+=temp_s/(2.0*gapcodons1*gapcodons2);
			 v_level[curr_degen]+=temp_v/(2.0*gapcodons1*gapcodons2);
			 
			 //cout<<" Degen: "<<curr_degen<<" Back S cnt: "<<temp_s<<" Back V cnt: "<<temp_v<<endl;
			 

			 break;
		 
		       case 2:
			 valid_steps=gapcodons1*gapcodons2*num_inters_wo_stops(codon1, codon2, num_diffs);
			 // cout<<"Codon 1: "<<codon1[0]<<codon1[1]<<codon1[2]<<" Codon 2: "<<codon2[0]<<codon2[1]<<codon2[2]
			 //  <<" Diffs: "<<num_diffs<<" Valid steps: "<<valid_steps;
			 for(l=0; l<2; l++)
			   {
			     intermediate_codon(codon1, codon2, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE)
			       {
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(codon1, diff1);
				 
				 curr_gen_code->count_by_degen(codon1, inter_codon, diff1, temp_s, temp_v); 
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;
				 
				 //cout<<" Inter codon: "<<inter_codon[0]<<inter_codon[1]<<inter_codon[2]<<" Degen: "<<curr_degen
				 // <<" C1->I S: "<<temp_s<<" C1->I V: "<<temp_v;

				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(inter_codon, diff2);
				 
				 curr_gen_code->count_by_degen(inter_codon, codon2, diff2, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;	
				 //cout<<" Degen: "<<curr_degen
				 //<<" I->C2 S: "<<temp_s<<" I->C2 V: "<<temp_v;
				  
			       }
		     
			     
			     intermediate_codon(codon2, codon1, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE)
			       {
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(codon2, diff1);			 
				 curr_gen_code->count_by_degen(codon2, inter_codon, diff1, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;							   
				 //cout<<" Inter codon(B): "<<inter_codon[0]<<inter_codon[1]<<inter_codon[2]<<" Degen: "<<curr_degen
				 //  <<"C2->I S: "<<temp_s<<" C2->I V: "<<temp_v;

				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(inter_codon, diff2);
				 curr_gen_code->count_by_degen(inter_codon, codon1, diff2, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps; 
				 //cout<<" Degen: "<<curr_degen
				   //  <<"I->C1 S: "<<temp_s<<" I->C1 V: "<<temp_v;
				
			       }		
			   }
			 //cout<<endl;
			 break;
			 
		       case 3:  
			 valid_steps=gapcodons1*gapcodons2*num_inters_wo_stops(codon1, codon2, num_diffs);
			 
			 for(l=0; l<6; l++)
			   {
			     intermediate_codon(codon1, codon2, num_diffs, l,  diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE &&
				 curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon2))==FALSE)
			       {
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(codon1, diff1);
				 curr_gen_code->count_by_degen(codon1, inter_codon, diff1, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;							   
				 
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(inter_codon, diff2);
				 curr_gen_code->count_by_degen(inter_codon, inter_codon2, diff2, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;							    
				 
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(inter_codon2, diff3);
				 curr_gen_code->count_by_degen(inter_codon2, codon2, diff3, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;							    
			       }
			     
			     intermediate_codon(codon2, codon1, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE &&
				 curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon2))==FALSE)
			       {
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(codon2, diff1);
				 curr_gen_code->count_by_degen(codon2, inter_codon, diff1, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;	
				 
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(inter_codon, diff2);
				 curr_gen_code->count_by_degen(inter_codon, inter_codon2, diff2, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;						    
				 
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(inter_codon2, diff3);
				 curr_gen_code->count_by_degen(inter_codon2, codon1, diff3, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;
			       }							   
			   }	
		       }
		   }
	       }
	   }
	}
    }
  l_level[0]=l_level[0]/2.0;
  l_level[2]=l_level[2]/2.0;
  l_level[4]=l_level[4]/2.0;

  //cout<<"Synonymous distance: v level (2,4): "<<v_level[2]<<","<<v_level[4]<<" s level (2,4: "<<s_level[2]<<","<<s_level[4]
    //  <<" l level (2,4): "<<l_level[2]<<","<<l_level[4]<<endl;

  dists[taxa_1].t_distance[taxa_2]=dists[taxa_2].t_distance[taxa_1]=
    (l_level[2]*(-0.5*log(1-2*(s_level[2]/l_level[2])-(v_level[2]/l_level[2]))+0.25*log(1-2*(v_level[2]/l_level[2])))+ 
     l_level[4]*(-0.5*log(1-2*(s_level[4]/l_level[4])-(v_level[4]/l_level[4]))+0.25*log(1-2*(v_level[4]/l_level[4]))))/
    (l_level[2]+l_level[4])-0.5*log(1-2*(v_level[4]/l_level[4]));
}




L_93_Non_Synon_dist::L_93_Non_Synon_dist(): Codon_dist()
{
}

L_93_Non_Synon_dist::L_93_Non_Synon_dist(Sequence_dataset *the_seq, Exchange *cexchange): Codon_dist(the_seq, cexchange)
{
}


void L_93_Non_Synon_dist::pairwise_dist(int taxa_1, int taxa_2)
{
  int i, j, k, l, codon1[3], codon2[3], inter_codon[3], inter_codon2[3], num_diffs, diff1, diff2, 
    diff3, curr_degen, numgap1, numgap2, gapcodons1, gapcodons2, gap1a, gap2a, gap1b, gap2b; 
  double l_level[5], s_level[5], v_level[5], temp_s, temp_v, valid_steps;


  for(i=0; i<4; i++)
    l_level[i]=s_level[i]=v_level[i]=0;

  for(i=0; i<((*curr_sequence)[0].Sequence_size()/3); i++)
    {
      gap_possible_codons(taxa_1, i, numgap1, gapcodons1, gap1a, gap2a);
      gap_possible_codons(taxa_2, i, numgap2, gapcodons2, gap1b, gap2b);
    

      for(j=0; j<3; j++)
	{
	  codon1[j]=(*curr_sequence)[taxa_1][3*i+j];
	  codon2[j]=(*curr_sequence)[taxa_2][3*i+j];
	}

      for (j=0; j<gapcodons1; j++)
	{
	  get_gap_codon(codon1, j, numgap1, gap1a, gap2a);

	  if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon1))==FALSE)
	   {      
	     for (k=0; k<gapcodons2; k++)
	       {
		 get_gap_codon(codon2, k, numgap2, gap1b, gap2b);
		 
		 if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon2))==FALSE)
		   {
		     num_diffs=curr_gen_code->multiple_subs(codon1, codon2);
		     
		     for(l=0; l<3; l++)
		       {
			 l_level[curr_gen_code->Li_93_degen(codon1, l)]+=1.0/(gapcodons1*gapcodons2);
			 l_level[curr_gen_code->Li_93_degen(codon2, l)]+=1.0/(gapcodons1*gapcodons2);
		       }
		     
		     switch (num_diffs)
		       {
		       case 1:
			 diff1=0;
			 while(codon1[diff1]==codon2[diff1])
			   diff1++;
			 
			 curr_degen=curr_gen_code->Li_93_degen(codon1, diff1);
			 
			 temp_s=temp_v=0;
			 curr_gen_code->count_by_degen(codon1, codon2, diff1, temp_s, temp_v);
			 
			 s_level[curr_degen]+=temp_s/(2.0*gapcodons1*gapcodons2);
			 v_level[curr_degen]+=temp_v/(2.0*gapcodons1*gapcodons2);
			 
			 curr_degen=curr_gen_code->Li_93_degen(codon2, diff1);
			 temp_s=temp_v=0;
			 
			 curr_gen_code->count_by_degen(codon2, codon1, diff1, temp_s, temp_v);
			 s_level[curr_degen]+=temp_s/(2.0*gapcodons1*gapcodons2);
			 v_level[curr_degen]+=temp_v/(2.0*gapcodons1*gapcodons2);
			 
			 break;
		 
		       case 2:
			 valid_steps=gapcodons1*gapcodons2*num_inters_wo_stops(codon1, codon2, num_diffs);
			 
			 for(l=0; l<2; l++)
			   {
			     intermediate_codon(codon1, codon2, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE)
			       {
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(codon1, diff1);
				 
				 curr_gen_code->count_by_degen(codon1, inter_codon, diff1, temp_s, temp_v); 
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;
				 
				 
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(inter_codon, diff2);
				 
				 curr_gen_code->count_by_degen(inter_codon, codon2, diff2, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;					  
			       }
		     
			     
			     intermediate_codon(codon2, codon1, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE)
			       {
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(codon2, diff1);
			 
				 curr_gen_code->count_by_degen(codon2, inter_codon, diff1, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;							   
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(inter_codon, diff2);
				 curr_gen_code->count_by_degen(inter_codon, codon1, diff2, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;
			       }		
			   }
			 break;
			 
		       case 3:  
			 valid_steps=gapcodons1*gapcodons2*num_inters_wo_stops(codon1, codon2, num_diffs);
			 
			 for(l=0; l<6; l++)
			   {
			     intermediate_codon(codon1, codon2, num_diffs, l,  diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE &&
				 curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon2))==FALSE)
			       {
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(codon1, diff1);
				 curr_gen_code->count_by_degen(codon1, inter_codon, diff1, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;							   
				 
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(inter_codon, diff2);
				 curr_gen_code->count_by_degen(inter_codon, inter_codon2, diff2, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;							    
				 
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(inter_codon2, diff3);
				 curr_gen_code->count_by_degen(inter_codon2, codon2, diff3, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;							    
			       }
			     
			     intermediate_codon(codon2, codon1, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE &&
				 curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon2))==FALSE)
			       {
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(codon2, diff1);
				 curr_gen_code->count_by_degen(codon2, inter_codon, diff1, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;	
				 
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(inter_codon, diff2);
				 curr_gen_code->count_by_degen(inter_codon, inter_codon2, diff2, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;						    
				 
				 temp_s=temp_v=0;
				 curr_degen=curr_gen_code->Li_93_degen(inter_codon2, diff3);
				 curr_gen_code->count_by_degen(inter_codon2, codon1, diff3, temp_s, temp_v);
				 s_level[curr_degen]+=temp_s/valid_steps;
				 v_level[curr_degen]+=temp_v/valid_steps;
			       }							   
			   }	
		       }
		   }
	       }
	   }
	}
    } 
 

  l_level[0]=l_level[0]/2.0;
  l_level[2]=l_level[2]/2.0;
  l_level[4]=l_level[4]/2.0;

  //cout<<"Non-Synonymous distance: v level (0,2): "<<v_level[0]<<","<<v_level[2]<<" s level (0,2): "<<s_level[0]<<","<<s_level[2]
    //  <<" l level (0,2): "<<l_level[0]<<","<<l_level[2]<<endl;
 

  dists[taxa_1].t_distance[taxa_2]=dists[taxa_2].t_distance[taxa_1]=
    -0.5*log(1-2*(s_level[0]/l_level[0])-(v_level[0]/l_level[0]))+0.25*log(1-2*(v_level[0]/l_level[0]))+
    (l_level[0]*-0.5*log(1-2*(v_level[0]/l_level[0]))+l_level[2]*-0.5*log(1-2*(v_level[2]/l_level[2])))/
    (l_level[2]+l_level[0]);
 
}




C_95_Synon_dist::C_95_Synon_dist() : Codon_dist()
{
}

C_95_Synon_dist::C_95_Synon_dist(Sequence_dataset *the_seq, Exchange *cexchange) : Codon_dist(the_seq, cexchange)
{
}


void C_95_Synon_dist::pairwise_dist(int taxa_1, int taxa_2)
{
  int i, j, k,l, codon1[3], codon2[3], inter_codon[3], inter_codon2[3], num_diffs, diff1, diff2, 
    diff3, curr_degen,  numgap1, numgap2, gapcodons1, gapcodons2, gap1a, gap2a, gap1b, gap2b; 
  double l_level[5], s_level[5], v_level[5],  Qs, Qa, A2s, A4, valid_steps;


  for(i=0; i<5; i++)
    l_level[i]=s_level[i]=v_level[i]=0;


  for(i=0; i<((*curr_sequence)[0].Sequence_size()/3); i++)
    {
      gap_possible_codons(taxa_1, i, numgap1, gapcodons1, gap1a, gap2a);
      gap_possible_codons(taxa_2, i, numgap2, gapcodons2, gap1b, gap2b);
    

      for(j=0; j<3; j++)
	{
	  codon1[j]=(*curr_sequence)[taxa_1][3*i+j];
	  codon2[j]=(*curr_sequence)[taxa_2][3*i+j];
	}

      for (j=0; j<gapcodons1; j++)
	{
	  get_gap_codon(codon1, j, numgap1, gap1a, gap2a);

	  if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon1))==FALSE)
	   {      
	     for (k=0; k<gapcodons2; k++)
	       {
		 get_gap_codon(codon2, k, numgap2, gap1b, gap2b);
		 
		 if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon2))==FALSE)
		   {
		     num_diffs=curr_gen_code->multiple_subs(codon1, codon2);
		     
		     for(l=0; l<3; l++)
		       {
			 l_level[curr_gen_code->Comeron_95_degen(codon1, l)]+=1.0/(gapcodons1*gapcodons2);
			 l_level[curr_gen_code->Comeron_95_degen(codon2, l)]+=1.0/(gapcodons1*gapcodons2);
		       }
		     
		     switch (num_diffs)
		       {
		       case 1:
			 diff1=0;
			 while(codon1[diff1]==codon2[diff1])
			   diff1++;
			
			 curr_gen_code->Comeron_levels(codon1, codon2, diff1, s_level, v_level, 
									2.0*gapcodons1*gapcodons2, curr_gen_code->Comeron_95_degen(codon1, diff1));
			
			 curr_gen_code->Comeron_levels(codon2, codon1, diff1, s_level, v_level, 
									2.0*gapcodons1*gapcodons2, curr_gen_code->Comeron_95_degen(codon2, diff1));
			 
			 break;
		 
		       case 2:
			 valid_steps=gapcodons1*gapcodons2*num_inters_wo_stops(codon1, codon2, num_diffs);
			 
			 for(l=0; l<2; l++)
			   {
			     intermediate_codon(codon1, codon2, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE)
			       {
				 curr_gen_code->Comeron_levels(codon1, inter_codon, diff1, s_level, v_level, 
									valid_steps, curr_gen_code->Comeron_95_degen(codon1, diff1));
				
				 curr_gen_code->Comeron_levels(inter_codon, codon2, diff2, s_level, v_level, 
									valid_steps, curr_gen_code->Comeron_95_degen(inter_codon, diff2)); 	  
			       }
		     
			     
			     intermediate_codon(codon2, codon1, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE)
			       {
				 curr_gen_code->Comeron_levels(codon2, inter_codon, diff1, s_level, v_level, 
										valid_steps, curr_gen_code->Comeron_95_degen(codon2, diff1));
						
				 curr_gen_code->Comeron_levels(inter_codon, codon1, diff2, s_level, v_level, 
										valid_steps, curr_gen_code->Comeron_95_degen(inter_codon, diff2));
			       }		
			   }
			 break;
			 
		       case 3:  
			 valid_steps=gapcodons1*gapcodons2*num_inters_wo_stops(codon1, codon2, num_diffs);
			 
			 for(l=0; l<6; l++)
			   {
			     intermediate_codon(codon1, codon2, num_diffs, l,  diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE &&
				 curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon2))==FALSE)
			       {
				 curr_gen_code->Comeron_levels(codon1, inter_codon, diff1, s_level, v_level, 
										valid_steps, curr_gen_code->Comeron_95_degen(codon1, diff1));
			       					   	
				 curr_gen_code->Comeron_levels(inter_codon, inter_codon2, diff2, s_level, 
										v_level, valid_steps, 
										curr_gen_code->Comeron_95_degen(inter_codon, diff2));
	
				 curr_gen_code->Comeron_levels(inter_codon2, codon2, diff3, s_level, 
										v_level, valid_steps, 
										curr_gen_code->Comeron_95_degen(inter_codon2, diff3));  
			       }
			     
			     intermediate_codon(codon2, codon1, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE &&
				 curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon2))==FALSE)
			       {
				 curr_gen_code->Comeron_levels(codon2, inter_codon, diff1, s_level, v_level, 
									valid_steps, curr_gen_code->Comeron_95_degen(codon2, diff1));
			       					   	
				 curr_gen_code->Comeron_levels(inter_codon, inter_codon2, diff2, s_level, 
										v_level, valid_steps,
										curr_gen_code->Comeron_95_degen(inter_codon, diff2));
	
				 curr_gen_code->Comeron_levels(inter_codon2, codon1, diff3, s_level, 
										v_level, valid_steps, 
										curr_gen_code->Comeron_95_degen(inter_codon2, diff3));  
			       }							   
			   }	
		       }
		   }
	       }
	   }
	}
    }
  l_level[0]=l_level[0]/2.0;
  l_level[1]=l_level[1]/2.0;
  l_level[2]=l_level[2]/2.0;
  l_level[4]=l_level[4]/2.0;



  Qs=(v_level[1]+v_level[4])/(l_level[1]+l_level[4]);
  Qa=(v_level[0]+v_level[2])/(l_level[0]+l_level[2]);

  A2s=-0.5*log(1-2*(s_level[2]/l_level[2])-Qa)                      +0.25*log(1-2*Qa);
  A4 =-0.5*log(1-2*(s_level[4]/l_level[4])-(v_level[4]/l_level[4])) +0.25*log(1-2*(v_level[4]/l_level[4]));
 
  dists[taxa_1].t_distance[taxa_2]=dists[taxa_2].t_distance[taxa_1]=
   ((l_level[2]*A2s+l_level[4]*A4)/(l_level[2]+l_level[4]))-0.5*log(1-2*Qs);
}




C_95_Non_Synon_dist::C_95_Non_Synon_dist(): Codon_dist()
{
}

C_95_Non_Synon_dist::C_95_Non_Synon_dist(Sequence_dataset *the_seq, Exchange *cexchange): Codon_dist(the_seq, cexchange)
{
}


void C_95_Non_Synon_dist::pairwise_dist(int taxa_1, int taxa_2)
{
  int i, j, k, l, codon1[3], codon2[3], inter_codon[3], inter_codon2[3], num_diffs, diff1, diff2, 
    diff3, curr_degen, numgap1, numgap2, gapcodons1, gapcodons2, gap1a, gap2a, gap1b, gap2b; 
  double l_level[5], s_level[5], v_level[5], Qs, Qa, A0, A2v, valid_steps;


  for(i=0; i<4; i++)
    l_level[i]=s_level[i]=v_level[i]=0;



  for(i=0; i<((*curr_sequence)[0].Sequence_size()/3); i++)
    {
      gap_possible_codons(taxa_1, i, numgap1, gapcodons1, gap1a, gap2a);
      gap_possible_codons(taxa_2, i, numgap2, gapcodons2, gap1b, gap2b);
    

      for(j=0; j<3; j++)
	{
	  codon1[j]=(*curr_sequence)[taxa_1][3*i+j];
	  codon2[j]=(*curr_sequence)[taxa_2][3*i+j];
	}

      for (j=0; j<gapcodons1; j++)
	{
	  get_gap_codon(codon1, j, numgap1, gap1a, gap2a);

	  if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon1))==FALSE)
	   {      
	     for (k=0; k<gapcodons2; k++)
	       {
		 get_gap_codon(codon2, k, numgap2, gap1b, gap2b);
		 
		 if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(codon2))==FALSE)
		   {
		     num_diffs=curr_gen_code->multiple_subs(codon1, codon2);
		     
		     for(l=0; l<3; l++)
		       {
			 l_level[curr_gen_code->Comeron_95_degen(codon1, l)]+=1.0/(gapcodons1*gapcodons2);
			 l_level[curr_gen_code->Comeron_95_degen(codon2, l)]+=1.0/(gapcodons1*gapcodons2);
		       }
		     
		     switch (num_diffs)
		       {
		       case 1:
			 diff1=0;
			 while(codon1[diff1]==codon2[diff1])
			   diff1++;
			
			 curr_gen_code->Comeron_levels(codon1, codon2, diff1, s_level, v_level, 
									2.0*gapcodons1*gapcodons2, curr_gen_code->Comeron_95_degen(codon1, diff1));
			
			 curr_gen_code->Comeron_levels(codon2, codon1, diff1, s_level, v_level, 
									2.0*gapcodons1*gapcodons2, curr_gen_code->Comeron_95_degen(codon2, diff1));
			 
			 break;
		 
		       case 2:
			 valid_steps=gapcodons1*gapcodons2*num_inters_wo_stops(codon1, codon2, num_diffs);
			 
			 for(l=0; l<2; l++)
			   {
			     intermediate_codon(codon1, codon2, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE)
			       {
				 curr_gen_code->Comeron_levels(codon1, inter_codon, diff1, s_level, v_level, 
									valid_steps, curr_gen_code->Comeron_95_degen(codon1, diff1));
				
				 curr_gen_code->Comeron_levels(inter_codon, codon2, diff2, s_level, v_level, 
									valid_steps, curr_gen_code->Comeron_95_degen(codon1, diff2)); 	  
			       }
		     
			     
			     intermediate_codon(codon2, codon1, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE)
			       {
				 curr_gen_code->Comeron_levels(codon2, inter_codon, diff1, s_level, v_level, 
										valid_steps, curr_gen_code->Comeron_95_degen(codon2, diff1));
						
				 curr_gen_code->Comeron_levels(inter_codon, codon1, diff2, s_level, v_level, 
										valid_steps, curr_gen_code->Comeron_95_degen(inter_codon, diff2));
			       }		
			   }
			 break;
			 
		       case 3:  
			 valid_steps=gapcodons1*gapcodons2*num_inters_wo_stops(codon1, codon2, num_diffs);
			 
			 for(l=0; l<6; l++)
			   {
			     intermediate_codon(codon1, codon2, num_diffs, l,  diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE &&
				 curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon2))==FALSE)
			       {
				 curr_gen_code->Comeron_levels(codon1, inter_codon, diff1, s_level, v_level, 
										valid_steps, curr_gen_code->Comeron_95_degen(codon1, diff1));
			       					   	
				 curr_gen_code->Comeron_levels(inter_codon, inter_codon2, diff2, s_level, 
										v_level, valid_steps, 
										curr_gen_code->Comeron_95_degen(inter_codon, diff2));
	
				 curr_gen_code->Comeron_levels(inter_codon2, codon2, diff3, s_level, 
										v_level, valid_steps, 
										curr_gen_code->Comeron_95_degen(inter_codon2, diff3));  
			       }
			     
			     intermediate_codon(codon2, codon1, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE &&
				 curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon2))==FALSE)
			       {
				 curr_gen_code->Comeron_levels(codon2, inter_codon, diff1, s_level, v_level, 
									valid_steps, curr_gen_code->Comeron_95_degen(codon2, diff1));
			       					   	
				 curr_gen_code->Comeron_levels(inter_codon, inter_codon2, diff2, s_level, 
										v_level, valid_steps,
										curr_gen_code->Comeron_95_degen(inter_codon, diff2));
	
				 curr_gen_code->Comeron_levels(inter_codon2, codon1, diff3, s_level, 
										v_level, valid_steps, 
										curr_gen_code->Comeron_95_degen(inter_codon2, diff3));  
			       }							   
			   }	
		       }
#if 0
		     switch (num_diffs)
		       {
		       case 1:
			 diff1=0;
			 while(codon1[diff1]==codon2[diff1])
			   diff1++;
			 
			 curr_degen=curr_gen_code->Comeron_95_degen(codon1, diff1);
			 curr_gen_code->Comeron_levels(codon1, codon2, diff1, s_level[curr_degen], v_level[curr_degen], 
									2.0*gapcodons1*gapcodons2);
			 
			 curr_degen=curr_gen_code->Comeron_95_degen(codon2, diff1);
			 curr_gen_code->Comeron_levels(codon2, codon1, diff1, s_level[curr_degen], v_level[curr_degen], 
									2.0*gapcodons1*gapcodons2);			 
			 break;
		 
		       case 2:
			 valid_steps=gapcodons1*gapcodons2*num_inters_wo_stops(codon1, codon2, num_diffs);
			 
			 for(l=0; l<2; l++)
			   {
			     intermediate_codon(codon1, codon2, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE)
			       {
				 curr_degen=curr_gen_code->Comeron_95_degen(codon1, diff1);
				
				 curr_gen_code->Comeron_levels(codon1, inter_codon, diff1, s_level[curr_degen], 
										    v_level[curr_degen], valid_steps);
				 
				 curr_degen=curr_gen_code->Comeron_95_degen(inter_codon, diff2);				
				 curr_gen_code->Comeron_levels(inter_codon, codon2, diff2, s_level[curr_degen], 
										    v_level[curr_degen], valid_steps);			  
			       }
		     
			     
			     intermediate_codon(codon2, codon1, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE)
			       {
				 curr_degen=curr_gen_code->Comeron_95_degen(codon2, diff1);				
				 curr_gen_code->Comeron_levels(codon2, inter_codon, diff1, s_level[curr_degen], 
										    v_level[curr_degen], valid_steps);
				 
				 curr_degen=curr_gen_code->Comeron_95_degen(inter_codon, diff2);				
				 curr_gen_code->Comeron_levels(inter_codon, codon1, diff3, s_level[curr_degen], 
										    v_level[curr_degen], valid_steps);		
				
			       }		
			   }
			 break;
			 
		       case 3:  
			 valid_steps=gapcodons1*gapcodons2*num_inters_wo_stops(codon1, codon2, num_diffs);
			 
			 for(l=0; l<6; l++)
			   {
			     intermediate_codon(codon1, codon2, num_diffs, l,  diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE &&
				 curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon2))==FALSE)
			       {
				 curr_degen=curr_gen_code->Comeron_95_degen(codon1, diff1);		
				 curr_gen_code->Comeron_levels(codon1, inter_codon, diff1, s_level[curr_degen], v_level[curr_degen], 
									valid_steps);
			       					   
				 curr_degen=curr_gen_code->Comeron_95_degen(inter_codon, diff2);		
				 curr_gen_code->Comeron_levels(inter_codon, inter_codon2, diff2, s_level[curr_degen], 
										v_level[curr_degen], valid_steps);

				 curr_degen=curr_gen_code->Comeron_95_degen(inter_codon2, diff3);		
				 curr_gen_code->Comeron_levels(inter_codon2, codon2, diff3, s_level[curr_degen], 
										v_level[curr_degen], valid_steps);  
			       }
			     
			     intermediate_codon(codon2, codon1, num_diffs, l, diff1, diff2, diff3, inter_codon, inter_codon2);
			     
			     if (curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon))==FALSE &&
				 curr_gen_code->is_stop(curr_gen_code->get_codon_num(inter_codon2))==FALSE)
			       {
				 curr_degen=curr_gen_code->Comeron_95_degen(codon2, diff1);		
				 curr_gen_code->Comeron_levels(codon2, inter_codon, diff1, s_level[curr_degen], v_level[curr_degen], 
									valid_steps);
			       					   
				 curr_degen=curr_gen_code->Comeron_95_degen(inter_codon, diff2);		
				 curr_gen_code->Comeron_levels(inter_codon, inter_codon2, diff2, s_level[curr_degen], 
										v_level[curr_degen], valid_steps);

				 curr_degen=curr_gen_code->Comeron_95_degen(inter_codon2, diff3);		
				 curr_gen_code->Comeron_levels(inter_codon2, codon1, diff3, s_level[curr_degen], 
										v_level[curr_degen], valid_steps);  
			       }							   
			   }	
		       }
#endif
		   }
	       }
	   }
	}
    }

 

  l_level[0]=l_level[0]/2.0;
  l_level[1]=l_level[1]/2.0;
  l_level[2]=l_level[2]/2.0;
  l_level[4]=l_level[4]/2.0;

 

  Qs=(v_level[1]+v_level[4])/(l_level[1]+l_level[4]);
  Qa=(v_level[0]+v_level[2])/(l_level[0]+l_level[2]);

  A2v=-0.5*log(1-2*(s_level[1]/l_level[1])-Qs)                      +0.25*log(1-2*Qs);
  A0 =-0.5*log(1-2*(s_level[0]/l_level[0])-(v_level[0]/l_level[0])) +0.25*log(1-2*(v_level[0]/l_level[0]));
 
 dists[taxa_1].t_distance[taxa_2]=dists[taxa_2].t_distance[taxa_1]=
   ((l_level[1]*A2v+l_level[0]*A0)/(l_level[1]+l_level[0]))-0.5*log(1-2*Qa);

 
}





AA_Percent_dist::~AA_Percent_dist()
{
  int i;
  
    if (num_diffs !=0) {
        for (i=0; i<curr_sequence->Num_sequences(); i++)
            delete[] num_diffs[i];
        delete[] num_diffs;
    }
}


void AA_Percent_dist::pairwise_dist(int taxa_1, int taxa_2)
{
 int i;
   for(i=0; i<(*curr_sequence)[0].Sequence_size(); i++)
     {
       if((*curr_sequence)[taxa_1][i] != 
	  (*curr_sequence)[taxa_2][i])
	 {
	   num_diffs[taxa_1][taxa_2]++;
	   num_diffs[taxa_2][taxa_1]++;
	 }
       
     }

   dists[taxa_1].t_distance[taxa_2]=(1.0*num_diffs[taxa_1][taxa_2])/
     (*curr_sequence)[0].Sequence_size();
   dists[taxa_2].t_distance[taxa_1]=(1.0*num_diffs[taxa_1][taxa_2])/
     (*curr_sequence)[0].Sequence_size();
}

