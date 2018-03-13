// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *           Aflow COREY OSES - Duke University 2013-2017                  *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses
// corey.oses@duke.edu

#ifndef _AUROSTD_XCOMBOS_H_
#define _AUROSTD_XCOMBOS_H_

namespace aurostd {
  class xcombos {
    private:
      //NECESSARY PRIVATE CLASS METHODS - START
      void free();
      void copy(const xcombos& b);
      //NECESSARY END CLASS METHODS - END

      bool m_initialized;
      vector<int> m_input;
      int n_choices, m_choose;
      bool m_permutation; //otherwise combination
      bool m_sort;
      
      bool m_started;
      bool m_exhausted; //all possibilities explored

      //for permutations
      vector<int> m_current;

      //for combinations
      vector<int> m_p;
      int m_x,m_y;
      
      void incrementPermutation();
      void incrementCombinations();
      void initializeCombinationsP();
      void setCombinationsIncrementParameters();
      void initialize();
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      xcombos();
      xcombos(const vector<int>& vec,bool sort=TRUE);      //permutations
      xcombos(int choice_count,int choose_count); //combinations, m choose n, m=choice_count, n=choose_count
      xcombos(const xcombos& b);
      //constructors - END
      ~xcombos();
      const xcombos& operator=(const xcombos& other) ;
      xcombos& operator++();
      void reset(); //reset with same parameters
      void reset(vector<int> vec,bool sort=TRUE); //reset with permutations
      void reset(int choice_count,int choose_count); //reset with combinations
      vector<int> getCombo() const; //grab current possibility
      int getN() const; //grab n (total count)
      int getM() const; //grab m (choose)
      vector<uint> getIndices() const; //get which indicies are 1
      bool increment(); //nice wrapper for ++, returns m_exhausted
    //NECESSARY PUBLIC CLASS METHODS - STOP
  };
} // namespace aurostd

#endif  // _AUROSTD_XCOMBOS_H_
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2018           *
// *           Aflow COREY OSES - Duke University 2013-2017                  *
// *                                                                         *
// ***************************************************************************
