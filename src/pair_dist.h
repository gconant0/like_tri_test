#ifndef ___PAIR_DIST_H___
#define ___PAIR_DIST_H___

#include "gen_dna_funcs.h"
#include "read_seq.h"
#include "gen_code.h"
#include "exchange.h"

struct distances 
{
  int id, num_com;
  double *t_distance;
};


class Pair_dist {
 public:
  Pair_dist();
  Pair_dist(Sequence_dataset *the_seq, Exchange *cexchange);
	void reset_seqs(Sequence_dataset *new_seq);
  void find_dist ();
  double return_dist(int taxa_1, int taxa_2);
  double return_closest_pair(int &taxa_1, int &taxa_2, 
			     BOOL &joined_1, BOOL &joined_2);
  void collapse_matrix (int on_1, int on_2);
  void print_dists();
	void ignore_ambigs()    {ambigs_ignored=TRUE;};
  virtual ~Pair_dist();

 protected:
	int **num_diffs, total_remaining, *joined, pos;
	distances *dists;
	BOOL ambigs_ignored;
	Sequence_dataset *curr_sequence;
	Exchange *curr_exchange;
  
	BOOL is_joined(int taxa);
	virtual void handle_gaps(int taxa, int character, int &start, int &end)=0;
	virtual void allocate_diff_matrix()=0;
	virtual void pairwise_dist(int taxa_1, int taxa_2)=0;

};




class Nucleotide_dist : public Pair_dist
{
 public:
  Nucleotide_dist() :Pair_dist() {};
  Nucleotide_dist(Sequence_dataset *the_seq, Exchange *cexchange) :
    Pair_dist(the_seq, cexchange)
    {allocate_diff_matrix();};

	~Nucleotide_dist();

 protected:
  void handle_gaps(int taxa, int character, int &start, int &end);
  void allocate_diff_matrix();
};



class Codon_dist : public Pair_dist
{
 public:
  Codon_dist() {};
  Codon_dist(Sequence_dataset *the_seq, Exchange *cexchange);
  ~Codon_dist();

 protected:
  double possible_diffs, **total_diffs;
  Genetic_code *curr_gen_code;
	void handle_gaps(int taxa, int character, int &start, int &end) {};
  void gap_possible_codons(int taxa, int character, int &numgaps, 
			   int &num_gap_codons, int &gap1, int &gap2 );
  void get_gap_codon(int codon[3], int gcodonnum, int numgaps, int gap1, int gap2);
  void allocate_diff_matrix();
  int num_inters_wo_stops(int codon1[3], int codon2[3], int num_diffs);
  void intermediate_codon(int codon1[3], int codon2[3], 
			  int num_diffs, int intermed_num, int &diff1,
			  int &diff2, int &diff3,
			  int intermed_codon[3], int intermed_codon1[3]);
  void proper_degen(double *correct_s, double *correct_v, double &v_0, double &v_2, 
		    double &v_4, double &s_0, double &s_2, double &s_4, int degen);
};


class Amino_Acid_dist : public Pair_dist
{
 public:
  Amino_Acid_dist() :
    Pair_dist() {num_diffs=0;};
  Amino_Acid_dist(Sequence_dataset *the_seq, Exchange *cexchange) :
    Pair_dist(the_seq, cexchange)
    {allocate_diff_matrix();};
 protected: 
  void handle_gaps(int taxa, int character, int &start, int &end);
  void allocate_diff_matrix();
};



class Percent_dist : public Nucleotide_dist
{
 public:
  Percent_dist() : Nucleotide_dist() {};
  Percent_dist(Sequence_dataset *the_seq, Exchange *cexchange) :
    Nucleotide_dist(the_seq, cexchange)
    {};
	~Percent_dist();
 protected:
  void pairwise_dist(int taxa_1, int taxa_2);
};



class JC_dist : public Nucleotide_dist 
{
 public:
  JC_dist();
  JC_dist(Sequence_dataset *the_seq, Exchange *cexchange);
 protected:
  void pairwise_dist(int taxa_1, int taxa_2);
};



class K2P_dist : public Nucleotide_dist 
{
 protected:
  void pairwise_dist(int taxa_1, int taxa_2);
};



class NG_86_Synon_dist : public Codon_dist 
{
 public:
  NG_86_Synon_dist();
  NG_86_Synon_dist(Sequence_dataset *the_seq, Exchange *cexchange);
 protected:
  void pairwise_dist(int taxa_1, int taxa_2);
};



class NG_86_Non_Synon_dist : public Codon_dist 
{
 public: 
  NG_86_Non_Synon_dist();
  NG_86_Non_Synon_dist(Sequence_dataset *the_seq, Exchange *cexchange);
 protected:
  void pairwise_dist(int taxa_1, int taxa_2);
};



class L_93_Synon_dist : public Codon_dist 
{
 public:
  L_93_Synon_dist();
  L_93_Synon_dist(Sequence_dataset *the_seq, Exchange *cexchange);
 protected:
  void pairwise_dist(int taxa_1, int taxa_2);
};



class L_93_Non_Synon_dist : public Codon_dist 
{
 public:
  L_93_Non_Synon_dist();
  L_93_Non_Synon_dist(Sequence_dataset *the_seq, Exchange *cexchange);
 protected:
  void pairwise_dist(int taxa_1, int taxa_2);
};

class C_95_Synon_dist : public Codon_dist 
{
 public:
  C_95_Synon_dist();
  C_95_Synon_dist(Sequence_dataset *the_seq, Exchange *cexchange);
 protected: 
  void pairwise_dist(int taxa_1, int taxa_2);
};



class C_95_Non_Synon_dist : public Codon_dist 
{
 public:
  C_95_Non_Synon_dist();
  C_95_Non_Synon_dist(Sequence_dataset *the_seq, Exchange *cexchange);
 protected:
  void pairwise_dist(int taxa_1, int taxa_2);
};

class AA_Percent_dist: public Amino_Acid_dist
{
 public:
  AA_Percent_dist() : Amino_Acid_dist() {};
  AA_Percent_dist(Sequence_dataset *the_seq, Exchange *cexchange):
    Amino_Acid_dist (the_seq, cexchange) 
    {};
  ~AA_Percent_dist();

protected:
  void pairwise_dist(int taxa_1, int taxa_2);
};

#endif



