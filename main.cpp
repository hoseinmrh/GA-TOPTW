#include "Interface.h"
#include <chrono>
using namespace std;
using namespace std::chrono;
using namespace std;
int main() {
    File file("rc208");
    file.read_file();
    float totalProfit = 0;
    totalProfit = calculate_profit();
    cout << "Total profit is --> " << totalProfit << '\n';
    for(int i = 0; i < 3; i++){
        int N = file.get_N();
        int V = file.get_V();
        auto start = high_resolution_clock::now();
        GA ga(N, V, 200,500,0.7);
        ga.genetic_algorithm();
//        ga.genetic_algorithm_time(60);
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<seconds>(stop - start);
        cout << "**************************************************" << '\n';
        cout << "Time taken: "
             << duration.count() << " seconds" << endl;

        cout<<"Iteration "<<i<<" finished!"<<'\n';
    }
    return 0;
}
