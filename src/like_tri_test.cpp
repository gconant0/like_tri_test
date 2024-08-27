#include <iostream>
#include <math.h>
#include "maxlike.h"
#include "codon_like.h"
#include "exchange.h"
#include "tree.h"
#include "read_tree.h"
#include "read_seq.h"
#include "gen_code.h"
#include "powell.h"
#include "pair_dist.h"
#include "gen_dna_funcs.h"
#include <string>

using namespace::std;
/*This program calculates the likelihood of a three-taxa tree under the molecular clock assumption.  The first two taxa in the input file are constrained to have the same branch lengths, while the third is the outgroup*/

void parse_args(int argc, char **argv, double &tol, BOOL &removegaps, BOOL &codon_freqs, BOOL &fixed_freqs, BOOL &mol_clock, CLOCK_TYPE &clock, int kaks_fix_branch[2], BOOL &branchkaks1, Exchange *curr_exchange, int &num_pns);



int main(int argc, char *argv[])
{
  int i=0, j,k, avg_cnt=0, ntaxa, nchars, temp_codon[3], kaks_fix_branch[2], num_pns;
  char datafile[100];
  double **p_ns, dist[3], final_lnL, temp, diff, codon_basefreqs[3][4], tol, old_b,  old_ka, old_ks, old_kaks, sat_lnL, unsat_lnL;;
  BOOL removegaps, codon_freqs, fixed_freqs, interleave, mol_clock, brnopt, all_ks_ok=TRUE, entered_initials=FALSE, branchkaks1=FALSE;
  CLOCK_TYPE clock;

  //Classes
  Exchange current_exchange;
  Three_Taxa_Tree get_tree;
  Tree *current_tree;
  Read_Sequence *read_seqs;
  Sequence_dataset *nuc_data, *new_data,  *current_data, *third_pos_data;
  Genetic_code *current_genetic_code;
  Like_model *model;
  Powell current_powell;
  L_93_Synon_dist *syn_dist;
  L_93_Non_Synon_dist *nonsyn_dist;
 
  if (argc>2)
  { 
    parse_args(argc, argv, tol, removegaps, codon_freqs, fixed_freqs, mol_clock, clock, kaks_fix_branch, branchkaks1, &current_exchange, num_pns);

    p_ns=new double*[num_pns];
    for(i=0; i<num_pns; i++)
      p_ns[i]=new double[1];

    if (mol_clock == TRUE) {
      current_exchange.set_mol_clock_3();
      current_exchange.set_clock_type(clock);
      
      brnopt=FALSE;
      
      cout<<"Molecular clock is ON\n";
      switch (clock) {
      case KS_CLOCK:
	cout<<"Ks values constrained for duplicates\n";
	break;
      case KA_CLOCK:
	cout<<"Ka values constrained for duplicates\n";
	break;
      }
     
    }
    else {
       brnopt=TRUE;
      
      cout<<"Molecular clock is OFF\n";
	
    }
    if (num_pns==2)
      cout<<"Ka/Ks values constrained for duplicates\n";
	
    if (removegaps ==TRUE)
      cout<<"Ignoring gap characters\n";
    else
      cout<<"Gap characters not ignored\n";
    
    if (argc>7) {
      cout<<"Reading params from command line: "<<argv[3][0]<<endl;
      j=3;
      while (argv[j][0] == '-') j++;

      if (argc<j+6) {
	cerr<<"Invalid number of arguments entered\n";
      }
      else { 
	entered_initials=TRUE;
	dist[0]=string_to_float(argv[j]);
	p_ns[0][0]=string_to_float(argv[j+1]);
	if (mol_clock==FALSE)
	  dist[1]=string_to_float(argv[j+2]);
	else
	  dist[1]=dist[0];
	p_ns[1][0]=string_to_float(argv[j+3]);
	dist[2]=string_to_float(argv[j+4]);
	p_ns[2][0]=string_to_float(argv[j+5]);
	
	cout<<"Using entered initial values: "<<dist[0]<<" "<<dist[1]<<" "<<dist[2]<<", "<<p_ns[0][0]<<" "<<p_ns[1][0]<<" "<<p_ns[2][0]<<endl;
      }
    }
    
    current_genetic_code=new Genetic_code(&current_exchange, TRUE, "\0");
    
    if (fixed_freqs==TRUE) 
      current_exchange.set_model(MG_94_K2P);                            //Muse and Gaut 94 model
    else
      current_exchange.set_model(MG_94_HKY);
    
    current_exchange.set_num_rates(1);
   
    //Gets the data

    switch(guess_dataformat(current_exchange.get_datafile(), strlen(current_exchange.get_datafile())))
       {
       case NEXUS:
         read_seqs=new Read_Nexus;
         break;
       case PIR:
         read_seqs=new Read_PIR;
         break;
       case PHYLIP: 
	 interleave=test_interleave(current_exchange.get_datafile());

	  if (interleave == TRUE)
	    read_seqs=new Read_Phylip_interleave;
	  else
	    read_seqs=new Read_Phylip_noninterleave;
	
         break;
       case FASTA:
         read_seqs=new Read_FASTA;
         break;
       } 

    nuc_data=read_seqs->get_dataset(ntaxa, nchars, current_exchange.get_datafile(), FALSE);
    if(removegaps == TRUE) {
      new_data= remove_gaps (nuc_data, FALSE);
      delete nuc_data;
      nuc_data=new_data;
      nchars=(*nuc_data)[0].Sequence_size();
      current_exchange.set_num_sites(nchars);
    }
       

    if (nuc_data==0)
      {
	delete read_seqs;
	return(-1);
      }

    current_exchange.set_num_taxa(ntaxa);
    current_exchange.set_num_sites(nchars);
    
    if (entered_initials ==FALSE ) {

      syn_dist= new L_93_Synon_dist(nuc_data, &current_exchange);
      nonsyn_dist= new L_93_Non_Synon_dist(nuc_data, &current_exchange);
      
      
      syn_dist->find_dist();
      nonsyn_dist->find_dist();
      
      
      
      dist[0]=dist[1]=syn_dist->return_dist(0,1);
      p_ns[0][0]=p_ns[1][0]=nonsyn_dist->return_dist(0,1)/syn_dist->return_dist(0,1);
      dist[2]=((syn_dist->return_dist(0,2)-(syn_dist->return_dist(0,1)/2))+(syn_dist->return_dist(1,2)-(syn_dist->return_dist(0,1)/2)))/2;
	  p_ns[num_pns-1][0]=((nonsyn_dist->return_dist(0,2)/syn_dist->return_dist(0,2))+(nonsyn_dist->return_dist(1,2)/syn_dist->return_dist(1,2)))/2;
	  
		
      
      
      
      
      
      if (!((dist[0] > 1e-300) &&(dist[0]<1e3))) {
	dist[0]=dist[1]=2.5;
      } 
      if (!((p_ns[0][0] > 1e-300) &&(p_ns[0][0]<1e20))) {
	p_ns[0][0]=p_ns[1][0]=0.1;
      }
      if (!((dist[2] > 1e-300) &&(dist[2]<1e20))) {
	dist[2]=2.5;
      } 
      if (!((p_ns[num_pns-1][0] > 1e-300) &&(p_ns[num_pns-1][0]<1e20))) {
	p_ns[num_pns-1][0]=0.1;
      }

    }
    if (fixed_freqs==FALSE)
      {
	//If we use non-equal basefreqs, we tell the optimizer not to optimize them
	current_exchange.fix_basefreq();

      
	if (codon_freqs==FALSE) {
	  for(i=0; i<4; i++) 
	    current_exchange.set_basefreqs(observed_basefreqs (nuc_data, i), i);  
	  //cout<<"Using base freqs: A: ";
	  //cout<<current_exchange.return_basefreq(0)<<" C: "<<current_exchange.return_basefreq(1)<<
	  //  " G: "<<current_exchange.return_basefreq(2)<<" T: "<<current_exchange.return_basefreq(3)<<endl;
	}
	else {
	  for(i=0; i<3; i++)
	    for(j=0; j<4; j++)
	      codon_basefreqs[i][j]= observed_codon_basefreqs(nuc_data, i, j);

	  current_exchange.set_use_codon_basefreqs(codon_basefreqs);
	  //cout<<"Using codon freqs: ";
	  // cout<<"Position 1: A: "<<current_exchange.return_codon_basefreq(0,0)<<" C: "<<current_exchange.return_codon_basefreq(0,1)<<
	  //  " G:"<<current_exchange.return_codon_basefreq(0,2)<<" T: "<<current_exchange.return_codon_basefreq(0,3)<<endl;
	  //cout<<"Position 2: A: "<<current_exchange.return_codon_basefreq(1,0)<<" C: "<<current_exchange.return_codon_basefreq(1,1)<<
	  //  " G: "<<current_exchange.return_codon_basefreq(1,2)<<" T: "<<current_exchange.return_codon_basefreq(1,3)<<endl;
	  //cout<<"Positon 3: A: "<<current_exchange.return_codon_basefreq(2,0)<<" C: "<<current_exchange.return_codon_basefreq(2,1)<<
	  //  " G: "<<current_exchange.return_codon_basefreq(2,2)<<" T: "<<current_exchange.return_codon_basefreq(2,3)<<endl;

	}
      }


    current_data = current_genetic_code->make_codon_sequence(nuc_data);
	
    
    //Gets the tree
    current_tree=get_tree.create_three_taxa_tree(&current_exchange, current_data, dist);
    
    if (current_tree==0)
      {
	delete read_seqs;
	delete current_data;
	delete syn_dist;
	delete nonsyn_dist;
	delete nuc_data;
	return(-1);
      }

    (*current_tree)[0]->set_p_nonsyn_num(0);
    if (num_pns==2) {
      (*current_tree)[1]->set_p_nonsyn_num(0);
      (*current_tree)[2]->set_p_nonsyn_num(1);
    }
    else {
      (*current_tree)[1]->set_p_nonsyn_num(1);
      (*current_tree)[2]->set_p_nonsyn_num(2);
    }
      
    current_exchange.set_num_p_nonsyn(num_pns, p_ns);
    
    if (branchkaks1 ==TRUE) {
      if (kaks_fix_branch[0] == 1) {
	cout<<"Fixing branch 1 to Ka/Ks = 1.0\n";
	current_exchange.set_branch_kaks1(0);
	current_exchange.set_p_non_syn(0, 1.0);
      }
      if (kaks_fix_branch[1] == 1) {
	cout<<"Fixing branch 2 to Ka/Ks = 1.0\n";
	current_exchange.set_branch_kaks1(1);
	current_exchange.set_p_non_syn(1, 1.0);
      }
      
    }
   
    switch(current_exchange.get_model())
      {
      case MG_94_K2P:	
		model = new MG_94_K2P_model (&current_exchange, current_data, current_tree, current_genetic_code);
		break;
      case MG_94_HKY:	
		model = new MG_94_HKY_model (&current_exchange, current_data, current_tree, current_genetic_code);
		break;
      }


   
    for(i=0; i<3; i++) {
      model->set_ratios((*current_tree)[i]);
      (*current_tree)[i]->set_brnlen(dist[i]/(*current_tree)[i]->get_syn_ratio());
    }

    if (mol_clock == TRUE) {
      if (clock == KS_CLOCK) {
	model->set_save_k_common(dist[2]);
      
	model->set_save_k_dupl(dist[0]);
      
	cout<<"Setting common: "<<dist[2]
	    <<" and split: "<<dist[0]<<endl;
      }
      else if (clock == KA_CLOCK) {
	cout<<"Initial: d:"<<dist[0]<<" brnlen: "<<(*current_tree)[0]->get_brnlen()<<" ratios: "<<(*current_tree)[0]->get_nonsyn_ratio()
	    <<", "<<(*current_tree)[0]->get_syn_ratio()<<" Pnonsyn: "<<current_exchange.get_p_non_syn((*current_tree)[0]->get_p_nonsyn_num())<<endl;
	cout<<"Initial: d:"<<dist[2]<<" brnlen: "<<(*current_tree)[2]->get_brnlen()<<" ratios: "<<(*current_tree)[2]->get_nonsyn_ratio()
	    <<", "<<(*current_tree)[2]->get_syn_ratio()<<" Pnonsyn: "<<current_exchange.get_p_non_syn((*current_tree)[1]->get_p_nonsyn_num())<<endl;
	model->set_save_k_common((*current_tree)[2]->get_brnlen()*(*current_tree)[2]->get_nonsyn_ratio());
	
	model->set_save_k_dupl((*current_tree)[0]->get_brnlen()*(*current_tree)[0]->get_nonsyn_ratio());
	
	cout<<"Ka fixed: Setting common: "<<dist[2]*p_ns[2][0]
	    <<" and split: "<<dist[0]*p_ns[0][0]<<endl;
      }
      cout<<"Reinitializing parameters\n";
      model->reinit_params();   
    }

    cout<<current_exchange.get_num_params()<<endl;
   
    if(tol!=0)
      current_powell.set_tolerance(tol);


    cout<<"Begining global optimization\n"<<flush;
    current_exchange.set_saved_lnL(current_powell.Init_min(model, &current_exchange, brnopt));
    
    model->set_expect_subs();
    
    cout<<"Final likelihood: "<<current_exchange.get_saved_lnL()<<endl;
    unsat_lnL= current_exchange.get_saved_lnL();
    
 
    cout<<"Name\tKs\tKa\tKa/Ks\n";
    for(i=0; i<current_exchange.get_num_branches(); i++)
      {
	cout<<(*current_tree)[i]->get_name()<<"\t"
	    <<(*current_tree)[i]->expect_subs_site()<<"\t"
	    <<(*current_tree)[i]->expect_nsyn_subs_site()<<"\t"
		<<current_exchange.get_p_non_syn((*current_tree)[i]->get_p_nonsyn_num())<<endl; 
      }


    for (i=0; i<3; i++) {
      old_b=(*current_tree)[i]->get_brnlen();   
      old_kaks=current_exchange.get_p_non_syn((*current_tree)[i]->get_p_nonsyn_num());

      old_ka=(*current_tree)[i]->expect_nsyn_subs_site();
      old_ks=(*current_tree)[i]->expect_subs_site();
      if ((old_ks > 1e-4) && (old_ka != 0.0)) {

		current_exchange.set_p_non_syn((*current_tree)[i]->get_p_nonsyn_num(), old_ka/(10*old_ks));				     
	
		model->set_ratios((*current_tree)[i]);
	
		(*current_tree)[i]->set_brnlen(10*old_ks/(*current_tree)[i]->get_syn_ratio());
	
	
		model->recalculate_transprobs();
		sat_lnL=model->find_ln_like_single_rate();
	
		//cout<<"Sat lnL: "<<sat_lnL<<endl;
		//cout<<"Unsat lnL: "<<unsat_lnL<<endl;
	
		if (sat_lnL>=unsat_lnL)
			cout<<"WARNING:  Synonymous substitutions on branch "<<i<<" are likely saturated.  Sat lnL: "<<sat_lnL<<" Unsat lnL: "<<unsat_lnL<<"\n";
	
		(*current_tree)[i]->set_brnlen(old_b);
		current_exchange.set_p_non_syn((*current_tree)[i]->get_p_nonsyn_num(), old_kaks);
      }
      
    }
     
   
   


    for(i=0; i<3; i++)
      delete[] p_ns[i];
    delete[] p_ns;
    //delete syn_dist;
    //delete nonsyn_dist;
    delete model;
    delete current_tree;
    delete nuc_data;
    delete current_data;
    delete read_seqs;
    delete current_genetic_code;

    return(0);
  } 
  else
    {
      cerr<<"Usage: like_tri_test <alignment file> -m:<OFF/Ks/Ka/KaKs> (-f:<1/2>) (-nogap) (-obsfreqs/-codonfreqs) (-t:<Tolerance>) (Ks1 Ka1 Ks2 Ka2 Ks3 Ka3)\n";
      return(-1);
    }
}  //End main




void parse_args(int argc, char **argv, double &tol, BOOL &removegaps, BOOL &codon_freqs, BOOL &fixed_freqs, BOOL &mol_clock, CLOCK_TYPE &clock, int kaks_fix_branch[2], BOOL &branchkaks1,  Exchange *curr_exchange, int &num_pns)
{
  int i,j, groups;
  char datafile[100], groups_file[100], tol_string[30];
  removegaps=FALSE;
  codon_freqs=FALSE;
  fixed_freqs=TRUE;
  curr_exchange->set_genetic_code(UNIVERSAL);
  clock=KS_CLOCK;

  kaks_fix_branch[0]=0;
  kaks_fix_branch[1]=0;

  tol=0;
  num_pns=3;
  mol_clock=FALSE;


  strcpy(datafile, argv[1]);
  curr_exchange->set_datafile(datafile);
  if ((argv[2][3] != 'o') && (argv[2][3] != 'O')) {
   
    if ((argv[2][4] == 'a') ||(argv[2][4] == 'A')) {
      if (strlen(argv[2])>5) {
	num_pns=2;
      }
    
      else  {
	mol_clock=TRUE;
	clock=KA_CLOCK;
      }
    }
    else {
      mol_clock=TRUE;
      clock=KS_CLOCK; }
	  
  
  }
  i=3;
  while ((i<argc) &&(argv[i][0] == '-'))
    {
      switch(argv[i][1]) {
      case 'n':
      case 'N':
	removegaps=TRUE;
	i++;
	break;
      case 'c':
      case 'C':
	fixed_freqs=FALSE;
	codon_freqs=TRUE;
	i++;
	break;
      case 'o':
      case 'O':
	fixed_freqs=FALSE;
	i++;
	break;
      case 'v':
      case 'V':
	switch (argv[i][3]) {
	case 'u':
	case 'U':
	  curr_exchange->set_genetic_code(UNIVERSAL);
	  break;
	case 'v':
	case 'V':
	  curr_exchange->set_genetic_code(VERT_MITO);
	  break;
	case 'y':
	case 'Y':
	  curr_exchange->set_genetic_code(YEAST_MITO);
	  break;
	case 'm':
	case 'M':
	  curr_exchange->set_genetic_code(MOLD_MITO);
	  break;
	case 'i':
	case 'I':
	  curr_exchange->set_genetic_code(INVERT_MITO);
	  break;
	case 'c':
	case 'C':
	  curr_exchange->set_genetic_code(CILIATE_NUC);
	  break;
	case 'e':
	case 'E':
	  curr_exchange->set_genetic_code(ECHINO_MITO);
	  break;
	}
	i++;
	break;
      case 'g':
      case 'G':
	curr_exchange->set_model(C_00_K2P);
	groups=string_to_int(argv[i]);
	strcpy(groups_file, argv[i+1]);
	
	curr_exchange->curr_groups=new Amino_acid_group(groups, groups_file);
	i+=2;
	break;
      case 't':
      case 'T':
	j=3;
	while(j<strlen(argv[i]))
	  tol_string[j-3]=argv[i][j++];
	tol=string_to_float(tol_string);
	i++;
	break;
      case 'f':
      case 'F':
	if (mol_clock == FALSE) {
	branchkaks1=TRUE;
	
	if (argv[i][3] == '1')
	  kaks_fix_branch[0]=1;

	if (argv[i][3] == '2')
	  kaks_fix_branch[1]=1;
	}
	else
	  cout<<"ERROR: Cannot fix Ka/Ks ratio with molecular clock on\n";
	i++;
	break;
      }
    }
}
