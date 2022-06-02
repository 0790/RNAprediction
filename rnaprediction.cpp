#include<fstream>
#include<iostream>
#include<stdlib.h>
#include<cstring>
#include<cstdio>
#include<chrono>
#include<string>
/**
 * @file 
 * @authors Deepti Kumar (2018B5A70790H), Eva Tiwari (2018B5A70816H), Shreya Srirampur (2018B4A70886H), Anirudh A (2018B4A70936)
 *  
 * @section DESCRIPTION
 * The following program takes a character array representing the RNA sequence as input and prints the number of bases pairs that are present in the secondary structure of the sequence.
 * Rules for forming bases pairs are as follows:
 * 1. Pairs of bases match up; Each base matches with 1 other base.
 * 2. Adenine always matches with Uracil.
 * 3. Cytosine always matches with Guanine.
 * 4. There are no kinks in the folded molecule.
 * 5. Structures are knot-free. 
 * 6. No sharp turns allowed
 */

using namespace std;
using namespace std::chrono;
int static opt[100][100] ;
int l ;
char seq[100] ;
/**
 * @brief Structure to store the pair of valid bases
 * 
 */
struct Print{
    char a;
    char b;
    int a1 ;
    int b1 ;
};
Print printvalues[100] ;
/**
 * @brief A function to check whether the given base pair is possible or not based on rule 2 and 3
 * 
 * @param a Base at location t-1 in loop  
 * @param b Base at location j-1 in loop
 * @return int value is 1 if the basis form a pair, else 0
 */
int check_pair(char a, char b){
    if((a == 'A' && b == 'U' ) || (a=='U' && b == 'A') || (a=='C' && b == 'G')|| (a=='G' && b == 'C')) return 1 ; // according to rule 2 and 3
    return 0 ;
}
/**
 * @brief Function to check if basis is already in a pair, according to rule 1
 * 
 * @param i index of the current basis
 * @return int value is 1 if already in a pair, else 0
 */
int check_exists(int i){
    for(int k = 0 ; k < l ; k++){
        if(i == printvalues[k].a1 || i == printvalues[k].b1)return 1;
    }
    return 0;

}
/**
 * @brief Recursive function to traverse through the sequence and check the total number of pairs of bases possible based on the rules
 * 
 * @param i Index of the start of the sequrnce under consideration 
 * @param j Index of the end of the sequence under consideration
 * @return int number of bases pairs possible within the sequence starting from i and ending at j
 */
int OPT(int i, int j){
    if(i >= j - 4) return  0 ; // set the number of pairs to 0 if i >= j -4 according to rule 6
    if(opt[i][j] != -1) return opt[i][j] ; // value has been calculated, skip extra computation 
    int unpaired ;
    if(opt[i][j-1] == -1)opt[i][j -1] = OPT(i,j-1) ; 
    unpaired = opt[i][j-1] ; // extract and store in unpaired to check if there are more bases pairs if the current i and j are not paired
    
    int paired = 0 ;
    for( int  t = i ; t <= j-5 ; t++){
        int check = check_pair(seq[t-1] , seq[j-1]); // check if it is a valid pair according to rule 2 and 3
        
        if(check == 1) { // if it is a pair
            // calculate the number of pairs in sub problems and extract the maximum among them
            if(opt[i][t-1] == -1)opt[i][t-1] = OPT(i,t-1) ;
            if(opt[t+1][j-1] == -1) opt[t+1][j-1] = OPT(t+1 , j-1) ;
            int cal = 1 + opt[i][t-1] + opt[t+1][j-1] ;
        
            paired = max(paired ,cal ) ;
            
        }
    }  
    opt[i][j] = max(paired , unpaired) ; // check if it maximized when seq[i] and seq[j] are unpaired or paired
    
    return opt[i][j] ;
}
/**
 * @brief Function to initialise the Memoization matrix and to remove white spaces if any from the input squence
 * 
 * @param s Input sequence with or without spaces  
 * @return int The length of the sequence without spaces
 */
int init(string s){
    int n = 0 ;
    for(int i = 0 ; s[i] != '\0'  ; i++){
        if(s[i] != ' '){ //checks if the character in the string is whitespace or not
            seq[n] = s[i] ; // if not a whitespace, it is added to the char array seq
            n++ ;
        }
    }
    seq[n] = '\0' ;
    n++ ;
    memset(opt,-1,sizeof(opt)) ; // setting the 2D array with -1, to check if the pair of characters is unvisited
    for(int i = 0 ; i <= n ; i++){
        for(int j = 0 ; j <= n ; j++){
            if(i >= j -4 )opt[i][j] = 0 ; //Initialize the matrix with 0 for all pairs with i >= j-4 according to rule 6
        }

    }
    return n;
  
}
/**
 * @brief Fucntion to get the pairs of basis from the calculated maximum
 * 
 * @param i Index of the start of the sequrnce under consideration
 * @param j Index of the end of the sequence under consideration
 * @param l The length of the structure printvalues with the valid base pairs 
 */
void traceback(int i, int j){
    
    if( j <= i) return ;
   
    if (opt[i][j] == opt[i][j-1]) // if the basis in j is unpaired
    {
        
        
        traceback(i, j-1) ;
    }
    else{
        // if j forms a pair with index towards the left
        
        for(int k = i ; k < j - 4 ; k++){
            
            int check = check_pair(seq[k-1], seq[j-1]) ;
            if( check == 1){ 
                
                if(k-1 <0){
                    if(opt[i][j] == opt[k+1][j-1] +1){
                        int checkK = check_exists(k) ;
                        int checkJ = check_exists(j) ;
                        if(checkK == 0 && checkJ == 0){
                            printvalues[l].a = seq[k-1] ;
                            printvalues[l].a1 = k ;
                            printvalues[l].b = seq[j-1] ;
                            printvalues[l].b1 = j ;
                            l++ ;
                            
                        }
                        traceback(k+1, j-1) ;
                        
                    }
                }
                else if(opt[i][j] == opt[i][k-1] + opt[k+1][j-1] +1){
                    int checkK = check_exists(k) ;
                    int checkJ = check_exists(j) ;
                    if(checkK == 0 && checkJ == 0){
                        printvalues[l].a = seq[k-1] ;
                        printvalues[l].a1 = k ;
                        printvalues[l].b = seq[j-1] ;
                        printvalues[l].b1 = j ;
                        l++;  
                        
                    }
                    traceback(i, k-1) ;
                    traceback(k+1 , j-1) ;
                    
                }
            }
            
        }
    }

}
/**
 * @brief main function of the program, takes the RNA sequence as input annd calls the init and OPT function to calculae the number of pairs
 * 
 * @return int returns 0 if the run is successful.
 */
int main(){
   
    int n = 0 ;
    string s ;
    // Input the RNA sequence. Newline is the delimiter
    cout << "\nEnter sequence: " ;
    getline(cin ,s) ;
    auto start = high_resolution_clock::now() ; // to calculate the running time of the program - start
    
    n = init(s) ;

    // Call OPT function to calculate number of bases pairs
    for(int k = 5 ; k < n ; k++ ){
        for(int i = 1 ; i <= n-k ; i++){
            int j = k+ i ;
            opt[i][j] = OPT(i,j) ;
            
        }
    }
     
    traceback(1,n) ;
    auto stop = high_resolution_clock::now() ; // to calculate the running time of the program - stop
    auto duration = duration_cast<microseconds>(stop - start) ;
    // Print number of bases pairs.
    cout << "\nNumber of base pairs in the secondary structure is: "  << opt[1][n] << endl ;
    cout << "\nThe pairs are as follows:" << endl ;
   
    for(int  a = 0 ; a < l ; a++){
        cout<<"(" << printvalues[a].a << " , " << printvalues[a].b << ") at indices " << printvalues[a].a1 << " " << printvalues[a].b1 << endl;
    }
    
    cout << "\nRun time = " << duration.count() << " ms" ;
    
    return 0;
}