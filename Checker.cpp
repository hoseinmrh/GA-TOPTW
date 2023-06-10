#include "Interface.h"
using namespace std;

int main(){
    File file("2-pr05");
    file.read_file();
    float totalProfit = 0;
    totalProfit = calculate_profit();
    vector<int> checkingVector = {0, 170, 62, 91, 166, 107, 32, 30, 79, 144, 218, 178, 0, -1, 0, 34, 179, 18, 132, 143, 173, 85, 72, 49, 225, 0};
    TOP top(0);
    top.calculate_solution_final(checkingVector,1);
    cout<<"Finish"<<'\n';
    return 0;
}