//#include <iostream>
//#include <string>	
#include "Prova.h"
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <vector>
using namespace std;

vector<int> globale(100);

void modificaVettore ( vector<int> &nomeVec)
{
   for (int i = 0; i < nomeVec.size(); i++)
   {
	   nomeVec.at(i) = i;
   }
   return;
}

int main()
{
	string name;
	double ciao;

	Prova object1("Enrico");
	Prova object2("Loredana");

	object1.displayMessage();
	object2.displayMessage();

	cout << "Inserisci un nome da tastiera non più lungo di 5 caratteri:\n" << endl ;
	//cin >> name	;
	//object1.setName(name)	;
	object1.displayMessage();

	ciao = 6.73837572849386;
	cout << setprecision(3) << ciao;

	const int vecSize = 10;
	int vec[vecSize];
	int seed;

	for (int j=0; j<vecSize; j++)
	{
	   vec[j] = 1 + (rand() % 6);
	}

	srand(time(0));
	seed = rand();

	for (int i=0; i<10; i++)
	{
	   srand(i*seed);
	   for (int j=0; j < vecSize; j++)
	   {
	   cout << 1 + (rand() % 6) << setw(5); 
	   }

	   cout << endl;
	}

	vector<int> prova(100);
	vector<int> prova2(100);
	prova[10] = 20;

	prova2 = prova;	
	for (int i=0; i < prova.size()-70; i++)
	{
	   cout << prova[i] << setw(10) << prova2[i] << setw(10) << globale.at(i) << endl;
	}
	cout << endl;
	cout << prova.size() << endl;
	cout << prova.at(10) << endl;
	cout << prova2.at(10) << endl;
	prova2.at(20) = 30;
	prova.at(10) = 50;
	cout << prova.at(10) << endl;
	cout << prova2.at(10) << endl;
	modificaVettore(prova);
	modificaVettore(globale);
	for (int i=0; i < prova.size()-70; i++)
	{
	   cout << prova[i] << setw(10) << prova2[i] << setw(10) << globale.at(i) << endl;
	}
}

