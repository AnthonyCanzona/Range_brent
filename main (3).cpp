#include <iostream>
#include <cmath>
#include "range.hpp"

#define PI 3.14159265358979323846

void benchmark(std::function< double(double) > f, double start, double stop, double period, int bestEval, double trueAnswer);
double calculatePoints(double range, double actualRange){return std::max(0.0,std::min(10.0,-5-std::log10(std::abs(range-actualRange)/actualRange)));}
double calculateCompetitionPoints(int evals, int bestEval){return std::max(0.0 , std::min(10.0,10.0-2*std::log2((double)evals/bestEval)) );}


double f( double x ) {
    static uint32_t count{ 0 };
    if (std::isinf(x)) return count;
    ++count;
    return std::sin( x ); //sin(x)
}
double f1( double x ) {
    static uint32_t count{ 0 };
    if (std::isinf(x)) return count;
    ++count;
    return std::sin(x) + std::sin(2*x) + std::sin(3*x) + std::sin(4*x); // sin(x) + sin(2x) + sin(3x) + sin(4x) 
}
double f2( double x ) {
    static uint32_t count{ 0 };
    if (std::isinf(x)) return count;
    ++count;
    return std::sin(x) + std::sin(2*x+1); // sin(x) + sin(2x+1)
}
double f3( double x ) {
    static uint32_t count{ 0 };
    if (std::isinf(x)) return count;
    ++count;
    return std::sin(2*PI*x) + std::sin(3*x) + std::cos(8*PI*x); // sin(2πx) + sin(3x) + cos(8πx)
}
double f4( double x ) {
    static uint32_t count{ 0 };
    if (std::isinf(x)) return count;
    ++count;
    return std::sin(x*x) + std::sin(4*x); // sin(x²) + sin(4x)
}

/*
--------------------------------------------- SCOREBOARD --------------------------------------
NAME           | FUNCTION | RANGE DIFFERENCE | FUNCTION EVALS | METHOD | ITERATIONS | EPSILON |
---------------+----------+------------------+----------------+--------+------------+---------+
               |          |                  |                |        |            |         |
Winston        | f        | 0.00000000e+00   | 37             | Brent  | 10         | 1e-14   |
Winston        | f1       | 0.00000000e+00   | 117            | Brent  | 10         | 1e-14   |
Winston        | f2       | 4.44089210e-16   | 57             | Brent  | 10         | 1e-14   |
Winston        | f3       | 8.88178420e-16   | 625            | Brent  | 11         | 1e-14   |
Winston        | f4       | 4.44089210e-16   | 464            | Brent  | 11         | 1e-14   |
Vikram         | f        | 0.00000000e+00   | 11             | SPI    |            |         |
Vikram         | f1       | 8.88178420e-16   | 105            | SPI    |            |         |
Vikram         | f2       | 4.44089210e-16   | 90             | SPI    |            |         |
Vikram         | f3       | 0.00000000e+00   | 235            | SPI    |            |         |
---------------+----------+------------------+----------------+--------+------------+---------+
** FOR ALL SUBMISSIONS, DO NOT MODIFY YOUR CODE BETWEEN SUBMISSIONS TO GET BETTER SCORES 
** Dont try to specialize in a single function, this is just to compare the general algorithm with others results

**Note, the range difference may vary for some functions since Matlab is also slightly inaccurate
**      The answers coded in should be accurate to about 11 decimal points
*/

int main(){
/*  FUNCTIONS:
    NAME    RANGE       FUNCTION                                Notes:
    f:      [0,10]      sin(x)                                  regular function
    f1:     [0,10]      sin(x) + sin(2x) + sin(3x) + sin(4x)
    f2:     [0,10]      sin(x) + sin(2x+1)
    f3:     [0,10]      sin(2πx) + sin(3x) + cos(8πx)           High frequency test
    f4:     [0,10]      sin(x²) + sin(4x)                       Range gets marginally larger as you approach 10
*/
    //           [a,b]   max freq       best score so far   true answer (dont touch unless you are adding your own)
    benchmark(f , 0,10 , 1.0/(2.0*PI) ,       11         ,  2.0);               
    benchmark(f1, 0,10 ,    2.0/PI    ,      105         ,  6.464769742222904);
    benchmark(f2, 0,10 ,    1.0/PI    ,       57         ,  3.356512126558125);
    benchmark(f3, 0,10 ,      4       ,      625         ,  5.406391588667142); //min:-2.721193923126904 at x=7.86804460
                                                                                //max: 2.685197665540238 at x=0.253244272
    //Pi as the frequency is just an estimate since its in terms of x²                                                                            
    benchmark(f4, 0,10 ,      PI      ,      464         ,  3.977505083289189); //Min: -1.983859975667271 at x=7.417849720
                                                                                //Max: 1.993645107621918 at x=9.789847300
    
    //                                     4.461676814830946
    //endpoint test, min on left (pi/6,0.5), max on right (pi/2, 1)
    // benchmark(f ,PI/6,PI/2, 1.0/(2.0*PI),     49         ,  0.5 ); 

    return 0;
}

void benchmark(std::function< double(double) > f, double start, double stop, double maxfreq, int bestEval, double trueAnswer){
        double ans = range(f, start ,stop , maxfreq);
        const int totalEvals = (int)f(INFINITY);
        printf("Range: %.15f. Difference: %.8e\n",ans,std::abs(ans - trueAnswer));
        printf("Accuracy points: %.2f / 10\n",calculatePoints(ans, trueAnswer));
        printf("Competition points (Lowest is %d evals): %.2f / 10\n",bestEval,calculateCompetitionPoints(totalEvals, bestEval));
        std::cout << "Function evals: " << totalEvals << std::endl << std::endl;
}