//
// Created by Hosein on 5/16/2023.
//

#ifndef HW2_C___INTERFACE_H
#define HW2_C___INTERFACE_H
#include <time.h>
#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <cstdlib>
#define MAX 20
using namespace std;

class Vertex;
class TOP;
static vector<Vertex> vertexVector;

void add_to_vertex_vector(Vertex v);

class Vertex{
protected:
    int i; //vertex number
    float x; //x coordinate
    float y; //y coordinate
    float d; //service duration
    float profit; //profit of the duration
    float opening_time;
    float closing_time;

    int nothing;

public:
    Vertex(int a){
        nothing = a;
    }
    Vertex(int i_v, float x_v, float y_v, float d_v, float profit_v, float openingTime_v, float closingTime_v){
        i = i_v;
        x = x_v;
        y = y_v;
        d = d_v;
        profit = profit_v;
        opening_time = openingTime_v;
        closing_time = closingTime_v;
    }

    int getI() const {
        return i;
    }

    float getX() const {
        return x;
    }

    float getY() const {
        return y;
    }

    float getD() const {
        return d;
    }

    float getProfit() const {
        return profit;
    }

    float getOpeningTime() const {
        return opening_time;
    }

    float getClosingTime() const {
        return closing_time;
    }

};

class File {
protected:
    int N, V;
    string fileName;
public:

    File(string file_name){
        fileName = file_name;
    }

    int get_N(){
        return N;
    }

    int get_V(){
        return V;
    }



    void set_N(int n){
        N = n;
    }

    void set_V(int v){
        V = v;
    }

    void read_file() {
        string firstPath = "C:\\HOSEIN\\jozve and tamrin\\8th Semester\\ADA\\HW2 ALL\\HW2 C++\\Instances\\";
        string format = ".txt";
        string filePath = firstPath + fileName + format;
        ifstream new_file;
        // Open a file to perform a write operation using a file object.
        cout << filePath << '\n';
        new_file.open(filePath);

        if (!new_file) {
            cout << "Can't open file!" << '\n';
            exit(1);
        }
        if (new_file.is_open()) {
            string sa;
            int lineCount = 1;
            // Read data from the file object and put it into a string.
            while (getline(new_file, sa)) {
                if (lineCount == 2) {
                    lineCount++;
                    continue;
                }
                split_lines(sa, lineCount);
                lineCount++;
            }
            new_file.close();

        }
    }

    void split_lines(string line, int lineCount){
        float array[MAX];
        fill_n(array,MAX,-1);
        if(lineCount == 1){

            stringstream lineStream(line);
            int index = 0;
            while (lineStream.good() && index < MAX)
            {

                lineStream >> array[index];
                index ++;
            }
            set_V(array[1]);
            set_N(array[2]);

        }
        else{
            stringstream lineStream(line);
            int index = 0;
            int row = lineCount - 3;
            while (lineStream.good() && index < MAX)
            {


                lineStream >> array[index];
                index ++;
            }
            create_vertex(array);

        }
    }

    void create_vertex(float array[]){
        int i = array[0];
        float x = array[1];
        float y = array[2];
        float d = array[3];
        float profit = array[4];
        int c_index = not_minus1_index(array);
        float opening_time = array[c_index - 1];
        float closing_time = array[c_index];
        Vertex vertex(i,x,y,d,profit,opening_time,closing_time);
        add_to_vertex_vector(vertex);

    }

    int not_minus1_index(float array[]) {
        for (int i = MAX - 1; i > -1; i--) {
            if (array[i] != -1)
                return (i);

        }
    }

};

class TOP{
protected:
    int N;
    int V;
    float time;
    float profit;
    vector<int> finalSolution;
    int nothing;

public:
    TOP(int a ){ //Temp Constructor
        finalSolution.clear();
        nothing = a;
        time = 0;
        profit = 0;
    }

    const vector<int> &getFinalSolution() const {
        return finalSolution;
    }

    TOP(int n, int v){
        N = n;
        V = v;
        time = 0;
        profit = 0;
    }

    float getProfit() const {
        return profit;
    }

    void addProfit(float profit){
        TOP::profit += profit;
    }

    void setTime(float time) {
        TOP::time = time;
    }

    void addTime(float time){
        TOP::time += time;
    }

    float distance_time(Vertex current, Vertex next){
        int x1 = current.getX();
        int y1 = current.getY();
        int x2 = next.getX();
        int y2 = next.getY();

        return sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2) * 1.0);
    }

    void calculate_solution(vector<int>solution , int check){
        profit = 0;
        float travel_Time = 0;
        int length = solution.size();
        int visitedNodes = 0;
        int index = 0;
        Vertex current = vertexVector[0];
        Vertex head = vertexVector[0];
        Vertex next = vertexVector[0];
        while (index < length){

            if(solution[index] < 0){
                setTime(0); // New path
                index++;
                continue;
            }
            else{

                current = vertexVector[solution[index]];
                int visitResult = canBeVisited(current);
                if(visitResult == 0){
                    head = current; // Current holder
                    cout<<"We visit "<<current.getI()<<'\n';
                    cout<<"Time Before visit "<<time<<'\n';
                    addTime(current.getD());
                    addProfit(current.getProfit());
                    cout<<"Time after visit "<<time<<'\n';
                    if(solution[index] != 0){
                        visitedNodes++;
                    }
                }
                else if(visitResult == -1){
                    cout<<"We can't visit "<<current.getI()<<'\n';
                    cout<<"Time in this node "<<time<<'\n';
                    addTime(-1 * travel_Time);
                }
                else{
                    head = current;
                    cout<<"We can visit "<<current.getI()<<"But after "<<visitResult<<'\n';
                    addTime(visitResult);
                    cout<<"Time Before visit "<<time<<'\n';
                    addTime(current.getD());
                    addProfit(current.getProfit());
                    cout<<"Time after visit "<<time<<'\n';
                    if(solution[index] != 0){
                        visitedNodes++;
                    }
                }

                if(index + 1 != length){
                    next = vertexVector[solution[index + 1]];
                    travel_Time = distance_time(head, next);
                    cout<<"Travel time from "<<head.getI()<<" to "<<next.getI()<<" "<<travel_Time<<'\n';
                    addTime(travel_Time);
                }
            }
            index++;
        }

        if(check == 1){
            cout<<"We visit "<<visitedNodes<<" nodes with profit "<<profit<<'\n';
        }

    }

    void calculate_solution2(vector<int>solution , int check){
        profit = 0;
        float travel_Time = 0;
        int length = solution.size();
        int visitedNodes = 0;
        int index = 0;
        Vertex current = vertexVector[0];
        Vertex head = vertexVector[0];
        Vertex next = vertexVector[0];
        while (index < length){
//            cout<<solution[index]<<'\n';
            bool condition1 = solution[index] < 0;
            bool condition2 = (-1 * solution[index]) %3 == 0;
//            cout<<"First condition"<< condition1<<'\n';
//            cout<<"Second Condition "<<condition2<<'\n';
            if(condition1 && condition2){
                setTime(0); // New path
                index++;
                continue;
            }
            else{
//                cout<<solution[index]<<'\n';
                if(solution[index] < 0){
                    solution[index] = 0;
                }
                current = vertexVector[solution[index]];
                int visitResult = canBeVisited(current);
                if(visitResult == 0){
                    head = current; // Current holder
//                    cout<<"We visit "<<current.getI()<<'\n';
//                    cout<<"Time Before visit "<<time<<'\n';
                    addTime(current.getD());
                    addProfit(current.getProfit());
//                    cout<<"Time after visit "<<time<<'\n';
                    if(solution[index] != 0){
                        visitedNodes++;
                    }
                }
                else if(visitResult == -1){
//                    cout<<"We can't visit "<<current.getI()<<'\n';
//                    cout<<"Time in this node "<<time<<'\n';
                    addTime(-1 * travel_Time);
                }
                else{
                    head = current;
//                    cout<<"We can visit "<<current.getI()<<"But after "<<visitResult<<'\n';
                    addTime(visitResult);
//                    cout<<"Time Before visit "<<time<<'\n';
                    addTime(current.getD());
                    addProfit(current.getProfit());
//                    cout<<"Time after visit "<<time<<'\n';
                    if(solution[index] != 0){
                        visitedNodes++;
                    }
                }

                if(index + 1 != length){
                    if(solution[index + 1] < 0){
                        int value_check = -1 * solution[index + 1];
                        if(value_check %3 != 0){
                            solution[index + 1] = 0;
                        }

                    }
                    next = vertexVector[solution[index + 1]];
                    travel_Time = distance_time(head, next);
//                    cout<<"Travel time from "<<head.getI()<<" to "<<next.getI()<<" "<<travel_Time<<'\n';
                    addTime(travel_Time);
                }
            }
//            print_vector(solution,"");
            index++;
        }

        if(check == 1){
            cout<<"We visit "<<visitedNodes<<" nodes with profit "<<profit<<'\n';
        }

    }

    void calculate_solution2_final(vector<int>solution , int check){
        profit = 0;
        float travel_Time = 0;
        int length = solution.size();
        int visitedNodes = 0;
        int index = 0;
        Vertex current = vertexVector[0];
        Vertex head = vertexVector[0];
        Vertex next = vertexVector[0];
        while (index < length){
            bool condition1 = solution[index] < 0;
            bool condition2 = (-1 * solution[index]) %3 == 0;
            if(condition1 && condition2){
                finalSolution.push_back(-1);
                setTime(0); // New path
                index++;
                continue;
            }
            else{
//                cout<<solution[index]<<'\n';
                if(solution[index] < 0){
                    solution[index] = 0;
                }
                current = vertexVector[solution[index]];
                int visitResult = canBeVisited(current);
                if(visitResult == 0){
                    head = current; // Current holder
//                    cout<<"We visit "<<current.getI()<<'\n';
//                    cout<<"Time Before visit "<<time<<'\n';
                    addTime(current.getD());
                    addProfit(current.getProfit());
//                    cout<<"Time after visit "<<time<<'\n';
                    finalSolution.push_back(current.getI());
                    if(solution[index] != 0){
                        visitedNodes++;
                    }
                }
                else if(visitResult == -1){
//                    cout<<"We can't visit "<<current.getI()<<'\n';
//                    cout<<"Time in this node "<<time<<'\n';
                    addTime(-1 * travel_Time);
                }
                else{
                    head = current;
//                    cout<<"We can visit "<<current.getI()<<"But after "<<visitResult<<'\n';
                    addTime(visitResult);
//                    cout<<"Time Before visit "<<time<<'\n';
                    addTime(current.getD());
                    addProfit(current.getProfit());
                    finalSolution.push_back(current.getI());
//                    cout<<"Time after visit "<<time<<'\n';
                    if(solution[index] != 0){
                        visitedNodes++;
                    }
                }

                if(index + 1 != length){
                    if(solution[index + 1] < 0){
                        int value_check = -1 * solution[index + 1];
                        if(value_check %3 != 0){
                            solution[index + 1] = 0;
                        }

                    }
                    next = vertexVector[solution[index + 1]];
                    travel_Time = distance_time(head, next);
//                    cout<<"Travel time from "<<head.getI()<<" to "<<next.getI()<<" "<<travel_Time<<'\n';
                    addTime(travel_Time);
                }
            }
//            print_vector(solution,"");
            index++;
        }

        if(check == 1){
            cout<<"We visit "<<visitedNodes<<" nodes with profit "<<profit<<'\n';
        }

    }


    void calculate_solution_final(vector<int>solution , int check){
        profit = 0;
        float travel_Time = 0;
        int length = solution.size();
        int visitedNodes = 0;
        int index = 0;
        Vertex current = vertexVector[0];
        Vertex head = vertexVector[0];
        Vertex next = vertexVector[0];
        while (index < length){

            if(solution[index] < 0){
                setTime(0); // New path
                finalSolution.push_back(solution[index]);
                index++;
                continue;
            }
            else{

                current = vertexVector[solution[index]];
                int visitResult = canBeVisited(current);
                if(visitResult == 0){
                    head = current; // Current holder
                    cout<<"We visit "<<current.getI()<<'\n';
                    cout<<"Time Before visit "<<time<<'\n';
                    addTime(current.getD());
                    addProfit(current.getProfit());
                    cout<<"Time after visit "<<time<<'\n';
                    finalSolution.push_back(current.getI());
                    if(solution[index] != 0){
                        visitedNodes++;
                    }
                }
                else if(visitResult == -1){
                    cout<<"We can't visit "<<current.getI()<<'\n';
                    cout<<"Time in this node "<<time<<'\n';
                    addTime(-1 * travel_Time);
                }
                else{
                    head = current;
                    cout<<"We can visit "<<current.getI()<<"But after "<<visitResult<<'\n';
                    addTime(visitResult);
                    cout<<"Time Before visit "<<time<<'\n';
                    addTime(current.getD());
                    addProfit(current.getProfit());
                    finalSolution.push_back(current.getI());
                    cout<<"Time after visit "<<time<<'\n';
                    if(solution[index] != 0){
                        visitedNodes++;
                    }
                }

                if(index + 1 != length){
                    next = vertexVector[solution[index + 1]];
                    travel_Time = distance_time(head, next);
                    cout<<"Travel time from "<<head.getI()<<" to "<<next.getI()<<" "<<travel_Time<<'\n';
                    addTime(travel_Time);
                }
            }
            index++;
        }

        if(check == 1){
            cout<<"We visit "<<visitedNodes<<" nodes with profit "<<profit<<'\n';
        }

    }

    int canBeVisited(Vertex v){
        if (time >= v.getOpeningTime() && time <= v.getClosingTime())
            return 0;
        else if (time < v.getOpeningTime())
            return ceill(v.getOpeningTime() - time);
        else{
            return -1;
        }
    }

    vector<int> random_solution_generator2(){
        vector<int> numbers;
        unsigned num = chrono::system_clock::now().time_since_epoch().count();
        for(int i = 1; i < N+1; i++)
            numbers.push_back(i);

        shuffle(numbers.begin(), numbers.end(),default_random_engine(num));

        return numbers;
    }

    vector<int> to_teams(vector<int> numbers){
        vector<int> solution;
        solution.push_back(-1); // We start from 0
        int neg_counter = -2;
        solution.push_back(numbers[0]);
        float nTest = N;
        float vTest = V;
        int division = ceil(nTest/vTest);
        for (int i = 1; i<N; i++){
            if ((i % division)==0){
                solution.push_back(neg_counter--); //Each path ends to 0
                solution.push_back(neg_counter--); // To separate path
                solution.push_back(neg_counter--); // next path start
            }
            solution.push_back(numbers[i]);
        }
        solution.push_back(neg_counter--); //Last part
        return solution;

    }
};

class GA{
protected:
    int N;
    int V;
    int population_size;
    int generations;//    float crossover_probability;
    float mutation_probability;
    vector<vector<int>> population_vector;
    vector<vector<int>> child_vector;
    vector<vector<int>> parents_children;

    vector<vector<int>> population_vector_teams;
    vector<vector<int>> child_vector_teams;
    vector<vector<int>> parents_children_teams;

    vector<int> solution;
public:
    GA(int n, int v, int populationSize, int generations, float mutationProbability) : N(n),
    V(v),
    population_size(populationSize),
    generations(generations),
    mutation_probability(mutationProbability) {}

    vector<int> swapVector(vector<int> vector, int index1, int index2){
        int tmp = vector[index1];
        vector[index1] = vector[index2];
        vector[index2] = tmp;
        return vector;
    }

    vector<int> insertionVector(vector<int> vector1, int index1 , int index2){
        vector<int> newVector = vector1;
        newVector.erase(newVector.begin() + index1);
        if(index2 > index1){
            auto position = newVector.begin() + index2 - 1;
            newVector.insert(position , vector1[index1]);
        }
        else{
            auto position = newVector.begin() + index2;
            newVector.insert(position , vector1[index1]);
        }
        return newVector;
    }

    bool localSearchSwap(){
        vector<int> tempSolution = solution;
        int solutionLength = tempSolution.size();
        float bestProfit = calculate2(solution , 0);
        int bestI;
        int bestJ;
        int flag = 0;
        for(int i = 0; i < solutionLength; i++){
            if(tempSolution[i] < 0){ //We don't want to swap these
                continue;
            }
            for(int j = i; j < solutionLength; j++){
                if(tempSolution[j] < 0){
                    continue;
                }
                tempSolution = swapVector(tempSolution, i, j);
                float localProfit = calculate2(tempSolution , 0);
                if(localProfit > bestProfit){
                    flag = 1;
                    bestProfit = localProfit;
                    bestI = i;
                    bestJ = j;
                }
                tempSolution.clear();
                tempSolution = solution;
            }
        }

        if (flag == 1) {
            solution = swapVector(solution, bestI, bestJ);

            return true;
        }

        return false;

    }

    bool localSearchInsertion(){
        vector<int> tempSolution = solution;
        int solutionLength = tempSolution.size();
        float bestProfit = calculate2(solution , 0);
        int bestI;
        int bestJ;
        int flag = 0;
        for(int i = 0; i < solutionLength; i++){
            if(tempSolution[i] < 0){ //We don't want to mess with these
                continue;
            }
            for(int j = i; j < solutionLength; j++){
                if(tempSolution[j] < 0){
                    continue;
                }
                tempSolution = insertionVector(tempSolution, i, j);
                float localProfit = calculate2(tempSolution , 0);
                if(localProfit > bestProfit){
                    flag = 1;
                    bestProfit = localProfit;
                    bestI = i;
                    bestJ = j;
                }
                tempSolution.clear();
                tempSolution = solution;
            }
        }

        if (flag == 1) {
            solution = insertionVector(solution, bestI, bestJ);
            return true;
        }

        return false;
    }

    void local_search(vector<int> best_I){
        solution.clear();
        solution = best_I;
        bool result_swap = localSearchSwap();
        bool result_insertion = localSearchInsertion();
        while(result_swap && result_insertion){
            result_swap = localSearchSwap();
            result_insertion = localSearchInsertion();
        }
    }

    static float calculate2(vector<int> someSolution , int check){
//        counter++;
        TOP top(0);
        top.calculate_solution2(someSolution, check);
        return top.getProfit();
    }

    float calculate2_final(vector<int> someSolution , int check){
        TOP top(0);
        top.calculate_solution2_final(someSolution, check);
        vector<int> final_solution = top.getFinalSolution();
        print_vector_final(final_solution, "Final");
        return top.getProfit();
    }

    void create_initial_solution(){
        population_vector.clear();
        child_vector.clear();
        parents_children.clear();
        TOP top(N,V);
        for(int i = 0; i < population_size; i++) {
            //Create random individual
            vector<int> individual = top.random_solution_generator2();
            vector<int> individual_to_team = top.to_teams(individual);
            population_vector.push_back(individual);
            population_vector_teams.push_back(individual_to_team);
        }

    }

    void create_population(vector<int> indices){
        TOP top(N,V);
        population_vector.clear();
        population_vector_teams.clear();
        int counter = 0;
        for(int &individual : indices){
            if(counter == population_size){
                break;
            }
            vector<int> individual_vector = parents_children[individual];
            population_vector.push_back(individual_vector);
            vector<int> individual_to_team = top.to_teams(individual_vector);
            population_vector_teams.push_back(individual_to_team);
            counter++;
        }
    }

    void print_vector(vector<int> v , string title){
        cout<<title<<'\n';
        for(int & element: v){
            cout<<element<<" ";
        }
        cout<<'\n';
    }

    void print_vector_final(vector<int> v , string title){
        cout<<title<<'\n';
        for(int & element: v){
            if(element == -1){
                cout<<'\n';
                continue;
            }
            cout<<element<<" ";
        }
        cout<<'\n';
    }

    void print_population(vector<vector<int>> V){
        int counter = 1;
        for(vector<int> p: V){
            cout<<"Individual "<<counter<<'\n';
            print_vector(p, "");
            calculate2(p,1);
            cout<<"******************************\n";
            counter++;
        }
    }

    static bool compare( const vector<int>&a , const vector<int>&b){
        return calculate2(a,0) > calculate2(b,0);
    }


    void sort_population(){
        sort(population_vector_teams.begin(), population_vector_teams.end(), compare);
    }

    vector<int> pmx_crossover(const vector<int>& parent1, const vector<int>& parent2) {
        int n = parent1.size();
        vector<int> offspring(n);
        // Choose two random positions for the crossover

        int start = rand() % n;
        int end = rand() % n;
        if (start > end) {
            swap(start, end);
        }

        // Copy the segment between the two positions from parent1 to offspring
        for (int i = start; i <= end; i++) {
            offspring[i] = parent1[i];
        }

        // Create a mapping between the values in the segment and their indices
        unordered_map<int, int> mapping;
        for (int i = start; i <= end; i++) {
            mapping[parent1[i]] = i;
        }

        // Fill the remaining positions in offspring using the corresponding values from parent2
        for (int i = 0; i < n; i++) {
            if (i < start || i > end) {
                int value = parent2[i];
                while (mapping.find(value) != mapping.end()) {
                    value = parent2[mapping[value]];
                }
                offspring[i] = value;
            }
        }

        return offspring;
    }

    vector<int> mutation_swap(vector<int> mySolution){
        vector<int> newSolution = mySolution;
        int randomIndex1 = vertexGenerator(mySolution);
        int randomIndex2 = vertexGenerator(mySolution);
//        cout<<"Swap"<<'\n';
//        cout<<mySolution[randomIndex1]<<"  "<<mySolution[randomIndex2]<<'\n';
        int tmp = newSolution[randomIndex1];
        newSolution[randomIndex1] = newSolution[randomIndex2];
        newSolution[randomIndex2] = tmp;
        return newSolution;
    }

    int vertexGenerator(vector<int> mySolution){

        int length = mySolution.size();
        int randomIndex = (rand() % length);
        int vertex = mySolution[randomIndex];
        while (vertex < 0){
            randomIndex = (rand() % length);
            vertex = mySolution[randomIndex];
        }
//        cout<<randomIndex<<"  "<<vertex<<'\n';
//        print_vector(mySolution, "");

        return randomIndex;

    }

    vector<int> mutation_insertion(vector<int> mySolution){

        vector<int> newSolution = mySolution;
        int randomIndex1 = vertexGenerator(mySolution);
        int randomIndex2 = vertexGenerator(mySolution);
//        cout<<"Insertion"<<'\n';
//        cout<<mySolution[randomIndex1]<<"  "<<mySolution[randomIndex2]<<'\n';
        newSolution.erase(newSolution.begin() + randomIndex1);
        if(randomIndex2 > randomIndex1){
            auto position = newSolution.begin() + randomIndex2 - 1;
            newSolution.insert(position , mySolution[randomIndex1]);
        }
        else{
            auto position = newSolution.begin() + randomIndex2;
            newSolution.insert(position , mySolution[randomIndex1]);
        }
        return newSolution;

    }

    vector<int> tournament_selection(const vector<vector<int>>& population, int tournament_size) {
        int n = population.size();
        vector<int> parent_indices;

        for (int i = 0; i < n; i++) {
            // Choose tournament_size random individuals from the population
            vector<int> tournament;
            for (int j = 0; j < tournament_size; j++) {
                int index = rand() % n;
                tournament.push_back(index);
            }

            // Find the individual with the highest fitness in the tournament
            int winner_index = tournament[0];
            int max_fitness = calculate2(population[tournament[0]],0); // Assumes a function called calculate_fitness that returns the fitness of an individual
            for (int j = 1; j < tournament_size; j++) {
                int index = tournament[j];
                int fitness = calculate2(population[index],0);
                if (fitness > max_fitness) {
                    winner_index = index;
                    max_fitness = fitness;
                }
            }

            // Add the index of the winner to the list of parent indices
            parent_indices.push_back(winner_index);
        }

        return parent_indices;
    }

    vector<int> roulette_wheel_selection(const vector<vector<int>>& population) {
        int n = population.size();
        vector<int> parent_indices;

        // Calculate the total fitness of the population
        double total_fitness = 0;
        for (auto individual : population) {
            double fitness = calculate2(individual,0); // Assumes a function called calculate_fitness that returns the fitness of an individual
            total_fitness += fitness;
        }

        // Calculate the relative fitness of each individual
        vector<double> relative_fitness(n);
        for (int i = 0; i < n; i++) {
            double fitness = calculate2(population[i],0);
            relative_fitness[i] = fitness / total_fitness;
        }

        // Calculate the cumulative distribution function
        vector<double> cdf(n);
        cdf[0] = relative_fitness[0];
        for (int i = 1; i < n; i++) {
            cdf[i] = cdf[i - 1] + relative_fitness[i];
        }

        // Select parents using the roulette wheel
        for (int i = 0; i < n; i++) {
            // Generate a random number between 0 and 1
            double r = (double) rand() / RAND_MAX;

            // Find the index of the individual whose slice of the roulette wheel contains the random number
            int index = -1;
            for (int j = 0; j < n; j++) {
                if (r <= cdf[j]) {
                    index = j;
                    break;
                }
            }

            // Add the index of the selected individual to the list of parent indices
            parent_indices.push_back(index);
        }

        return parent_indices;
    }

    void bearing(vector<int> parent_indices){
        TOP top(N,V);
        int length = parent_indices.size();
        int counter = 0;
        child_vector.clear();
        child_vector_teams.clear();
        for(int i = 0; i < length; i++){
            for (int j = 0; j < length; j++){
                if(counter == population_size * 2){
                    break;
                }
                if( parent_indices[i] == parent_indices[j] ){
                    continue;
                }
                vector<int> parent1 = population_vector[parent_indices[i]];
                vector<int> parent2 = population_vector[parent_indices[j]];
                vector<int> offspring = pmx_crossover(parent1, parent2);
                child_vector.push_back(offspring);
                vector<int> offspring_to_team = top.to_teams(offspring);
                child_vector_teams.push_back(offspring_to_team);
                counter++;
            }
        }
    }

    float random_0_1(){
        return rand() % 100;
    }

    void add_parent_child(){
        TOP top(N,V);
        parents_children.clear();
        parents_children_teams.clear();
        for(vector<int> p: population_vector){
            parents_children.push_back(p);
        }

        for(vector<int> ch: child_vector){
            parents_children.push_back(ch);
        }

        for(vector<int> p_t: population_vector_teams){
            parents_children_teams.push_back(p_t);
        }

        for(vector<int> ch_t: child_vector_teams){
            parents_children_teams.push_back(ch_t);
        }


    }

    void mutation(){
        for(int i = 0; i < parents_children.size(); i++){
            float random_value = random_0_1() / 100;
//            cout<<random_value<<" random value"<<'\n';
            if(random_value < mutation_probability){
//                cout<<"We mutate"<<'\n';
                // We can mutate index i
                for(int mutation_number = 0; mutation_number < 10; mutation_number++){
                    float second_random_value = random_0_1();
                    vector<int> individual = parents_children[i];
                    if(second_random_value < 50){
                        individual = mutation_swap(individual);
                    }
                    else{
                        individual = mutation_insertion(individual);
                    }
                    // Assign to parent children
                    parents_children[i] = individual;
                }

            }
        }


    }

    vector<int> final_selection(){

        //Define kind of selection
        float random_value = random_0_1();
        vector<int> indices;
        int length_all = parents_children_teams.size();
        if(random_value < 50){
            indices = tournament_selection(parents_children_teams, length_all);
        }
        else{
            indices = roulette_wheel_selection(parents_children_teams);
        }

        return indices;

    }

    void genetic_algorithm(){
        srand(time(0));
        create_initial_solution();
        for(int i = 0; i < generations; i++){
//            init();
            //Select parents
            float random_value = random_0_1();
            vector<int> parent_indices;
            parent_indices.clear();
            if(random_value < 50){
                parent_indices = tournament_selection(population_vector_teams,
                                                      (3 * population_size) / 4);
            }
            else{
                parent_indices = roulette_wheel_selection(population_vector_teams);
            }
            //Create new generation
            bearing(parent_indices);
            vector<int> indices;
            indices.clear();
            add_parent_child();
            mutation();
//            local_search();
            indices = final_selection();
            create_population(indices); //Generate new population

            if(i % 50 == 0){
                cout<<"Generation "<<i<<'\n';
//                cout<<parents_children_teams.size()<<" P CH"<<'\n';
//                cout<<child_vector_teams.size()<<" CH"<<'\n';
//                cout<<population_vector_teams.size()<<" P"<<'\n';
            }
        }

//       Give out best individual
        sort_population();
//        print_population(population_vector_teams);
        vector<int> best = population_vector_teams[0];
//        local_search(best);
//        print_vector(best, "Best");
//        calculate2(best,1);
        calculate2_final(best,1);
    }

    void genetic_algorithm_time(int seconds){
        srand(time(0));
        create_initial_solution();
        auto Start = std::chrono::high_resolution_clock::now();
        int counter = 0;
        while(1){
            auto End = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> Elapsed = End - Start;
            if (Elapsed.count() >= seconds * 1000)
                break;
            float random_value = random_0_1();
            vector<int> parent_indices;
            parent_indices.clear();
            if(random_value < 50){
                parent_indices = tournament_selection(population_vector_teams,
                                                      (3 * population_size) / 4);
            }
            else{
                parent_indices = roulette_wheel_selection(population_vector_teams);
            }
            //Create new generation
            bearing(parent_indices);
            vector<int> indices;
            indices.clear();
            add_parent_child();
            mutation();
//            local_search();
            indices = final_selection();
            create_population(indices); //Generate new population

//            if(i % 50 == 0){
//                cout<<"Generation "<<i<<'\n';
//                cout<<parents_children_teams.size()<<" P CH"<<'\n';
//                cout<<child_vector_teams.size()<<" CH"<<'\n';
//                cout<<population_vector_teams.size()<<" P"<<'\n';
//            }
        }
//       Give out best individual
        sort_population();
//        print_population(population_vector_teams);
        vector<int> best = population_vector_teams[0];
//        cout<<counter<<'\n';
//        local_search(best);
//        print_vector(best, "Best");
//        calculate2(best,1);
        calculate2_final(best,1);
    }

};

void add_to_vertex_vector(Vertex v){
    vertexVector.push_back(v);
}

void print_file_info(){
    for(Vertex & v: vertexVector){
        cout<<"I "<<v.getI()<<", X "<<v.getX()<<", Y "<<
        v.getY()<<", D "<<v.getD()<<", P "<<v.getProfit()<<", OT "<<v.getOpeningTime()<<", CT "<<v.getClosingTime()<<'\n';
    }
}

float calculate_profit(){
    int profitSum = 0;
    for(int i = 0 ; i< vertexVector.size(); i++){
        profitSum = profitSum + vertexVector[i].getProfit();
    }
    return profitSum;
}

#endif //HW2_C___INTERFACE_H
