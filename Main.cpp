#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <unordered_map>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <limits>





struct Term{
    double coefficient;
    double constant;

    Term(double coefficient, double constant){
        this->coefficient = coefficient;
        this->constant = constant;
    }

    double getConstant(){
        return constant;
    }

    void setConstant(double constant){
        this->constant = constant;
    }
};

struct Matrix
{
    private:
    //The underlying matrix data
    std::vector<std::vector<Term>> mat;

    //Number of Rows
    int rows;

    //Number of Colums
    int cols;

    bool isVectorContains(std::vector<int> vector, int value){
        for(int val : vector){
            if(val == value){
                return true;
            }
        }
        return false;
    }

    public:
    Matrix(std::vector<std::vector<Term>> mat){
        this->mat = mat;
        rows = mat.size();
        if(rows != 0){
            cols = mat.at(0).size();
        }
    }

    double calcDeterminant(){
        double determinant = 1;

        std::vector<int> pivotRows;
        if(rows != cols){
            throw std::runtime_error("The Matrix is not square");
        }

        //Find the pivot column
        for (int columnStart = 0; columnStart < cols; columnStart++){
            int leadingRow = -1;
            double pivotValue = 0;
            for(int row = 0; row < rows; row++){
                double value = mat.at(row).at(columnStart).getConstant();
                if(value != 0 && !isVectorContains(pivotRows, row)){
                    leadingRow = row;
                    pivotValue = value;
                    break;
                }
            }
            if(leadingRow != -1){
                pivotRows.push_back(leadingRow);

                std::cout << leadingRow << " Lead row " << std::endl;
                determinant *= pivotValue;
                //With the leading row determined we need to scale it to be 1
                for(int column = columnStart; column < cols; column++){
                    mat.at(leadingRow).at(column).setConstant(mat.at(leadingRow).at(column).getConstant()/pivotValue);
                }

                for(int row = 0; row < rows; row++){
                    if(row != leadingRow){
                        double leadingValue = mat.at(row).at(columnStart).getConstant();
                        std::cout << leadingValue << "leadingValue" << std::endl;
                        for(int column = columnStart; column < cols; column++){
                            mat.at(row).at(column).setConstant(mat.at(row).at(column).getConstant() - leadingValue * mat.at(leadingRow).at(column).getConstant());
                        }
                    }
                }
            }
            else{
                determinant *= 0;
            }
        }
        return determinant;
    }

    void printMatrix(){
        std::cout << "rows: " << rows << " cols: " << cols << std::endl; 
        for(std::vector<Term> row : mat){
            for(Term term : row){
                std::cout << term.getConstant() << ", ";
            }
            std::cout << std::endl;
        }
    }
};



class Polynomial{
    private:
    std::vector<double> polynomial;

    public:

    Polynomial(std::vector<double> polynomial){
        this->polynomial = polynomial;
    }

    Polynomial(){
        polynomial = std::vector<double>();
    }

    std::vector<double> getAllZeros(){
        Polynomial currentPolynomial = polynomial;

        std::vector<double> zeros;
        for(int i = 0; i < polynomial.size() - 1; i++){
            double zero = currentPolynomial.getZero();
            currentPolynomial = Polynomial(currentPolynomial.syntheticDivision(zero).getCoefficientVector());
            currentPolynomial.printPolynomial();
            std::cout << zero << " has a zero" << std::endl;
            zeros.push_back(zero);
        }

        return zeros;
    }

    Polynomial syntheticDivision(double divisor){
        std::vector<double> dividend;
        int value = polynomial.at(polynomial.size() - 1);
        std::cout << value << " intit vaue "<< std::endl;
        dividend.push_back(value);
        for(int i = (polynomial.size() - 2); i >= 1; i--){
            for(double d : dividend){
                std::cout << d << " dividene " << divisor << std::endl;
            }
            std::cout << value << " intit vaue " << divisor << "idk what the fuck is wrong " << value * divisor << std::endl;
            value = value * divisor;
            std::cout << value << " and  " << polynomial.at(i) << " i " << i << std::endl;

            value += polynomial.at(i);
            dividend.push_back(value);
        }
        for(double d : dividend){
            std::cout << d << " dividene" << std::endl;
        }
        /*
        for(double coeff : dividend){
            std::cout << coeff << std::endl;
        }
        */
        std::reverse(dividend.begin(), dividend.end());
        return Polynomial(dividend);
    }

    double getZero(){
        int numIterations = 100;
        double xValue = 0;
        for(int i = 0; i < numIterations; i++){
            xValue = xValue - (getValue(xValue)/getDerivative(xValue));
        }
        //std::cout << getValue(xValue) << " for zero " << xValue << std::endl;
        return xValue;
    }

    double getDerivative(double xValue){
        double step = 0.0000001;
        double deltaY = getValue(xValue + step) - getValue(xValue);    
        return deltaY/step;
    }

    double getValue(double xValue){
        int size = polynomial.size();
        double value = 0;
        for(int i = 0; i < size; i++){
            value += polynomial.at(i) * std::pow(xValue, i);
        }
        return value;
    }

    int getDegree(){
        return polynomial.size() - 1;
    }

    double getCoefficient(int degree){
        return polynomial.at(degree);
    }

    std::vector<double> getCoefficientVector(){
        return polynomial;
    }

    void setCoefficientVector(std::vector<double> coefficientVector){
        polynomial = coefficientVector;
    }

    void addPolynomial(Polynomial add){
        Polynomial large;
        Polynomial small;

        std::vector<double> sum;

        //Determine which Polynomial has the higher degree
        if(add.getDegree() > polynomial.size()-1){
            large = add;
            small = polynomial;
        }else{
            large = polynomial;
            small = add;
        }

        std::cout << "The large degree is " << large.getDegree() << std::endl;

        for(int degree = 0; degree <= large.getDegree(); degree++){
            int coeffSum = large.getCoefficient(degree);
            if(small.getDegree() >= degree){
                coeffSum += small.getCoefficient(degree);
            }
            sum.push_back(coeffSum);
        }

        polynomial = sum;
    }

    void multiplyPolynomial(Polynomial multiply){
        std::vector<Polynomial> subTerms;

        for(int scalarIndex = 0; scalarIndex < polynomial.size(); scalarIndex++){
            std::vector<double> tempTerm;
            for(int i = 0; i < scalarIndex; i++){
                tempTerm.push_back(0);//Do increase the degree of the polynomial accordingly
            }
            double scalar = polynomial.at(scalarIndex);
            for(double coefficient : multiply.getCoefficientVector()){
                tempTerm.push_back(scalar * coefficient);
            }
            subTerms.push_back(Polynomial(tempTerm));
        }

        Polynomial sumTerms(std::vector<double>{0});
        for(Polynomial subPoly : subTerms){
            sumTerms.addPolynomial(subPoly);
        }
        polynomial = sumTerms.getCoefficientVector();
    }

    void printPolynomial(){
        for(double coeff : polynomial){
            std::cout << coeff << std::endl;
        }
    }
};

class Rational{
    private:
    Polynomial numerator;

    Polynomial denominator;

    static bool isVectorContains(std::vector<double> vector, int value){
    for(int val : vector){
        if(val == value){
            return true;
        }
    }
    return false;
}


    public:
    Rational(Polynomial numerator){
        this->numerator = numerator;
        denominator = Polynomial{std::vector<double>{1}};
    }

    Rational(Polynomial numerator, Polynomial denominator){
        this->numerator = numerator;
        this->denominator = denominator;
    }

    void simplify(){
        std::vector<double> numeratorZeros = numerator.getAllZeros();
        std::vector<double> denominatorZeros = denominator.getAllZeros();
        for(double zero : numeratorZeros){
            std::cout << zero << " is a zero in the num" << std::endl;
        }


        for(double zero : denominatorZeros){
            std::cout << zero << " is a zero in the den" << std::endl;
            if(isVectorContains(numeratorZeros, zero)){
                numerator = numerator.syntheticDivision(zero);
                denominator = denominator.syntheticDivision(zero);
            }
        }
    }

    void printRational(){
        std::cout << "Numerator Polynomial" << std::endl;
        numerator.printPolynomial();
    
        std::cout << "Denominator Polynomial" << std::endl;
        denominator.printPolynomial();
    }

};


struct WebPage{
    private:
    //Search Term
    std::string keyword;

    //Name of the Website
    std::string websiteName;

    //Vector of the WebPages this WebPage Links to
    std::vector<WebPage> links;

    //Map of WebPage names to their index
    std::unordered_map<std::string, bool> webPageMap;

    public:
    //Constructor for the webpage
    WebPage(std::string websiteName, std::string keyword){
        this->websiteName = websiteName;
        this->keyword = keyword;
    }

    void setLinks(std::vector<WebPage> links){
        this->links = links;
        
        for(int i = 0; i < links.size(); i++){
            webPageMap[links.at(i).getWebsiteName()] = true;
        }
    }

    std::string getWebsiteName(){
        return websiteName;
    }



    std::string loadWebPage(){
        std::cout << "Welcome to " << websiteName << "!! Where would you like to go next" << std::endl;
        for(WebPage page : links){
            std::cout << "Type " << page.getWebsiteName() << " to go to " << page.getWebsiteName() << std::endl;
        }
        std::string nextPage;
        std::cin >> nextPage;
        return nextPage;
    }

    bool isWebsiteAvailable(std::string name){
        if(webPageMap.at(name) == true){
            return true;
        }
        return false;
    }
};

class SearchEngine{
    private:
    //vector of the webpages inside this search engine
    std::vector<WebPage>* pages;

    //Current WebPage
    int websiteID;

    //Map of WebPage names to their index
    std::unordered_map<std::string, int> webPageMap;

    public:
    //Constructor for the SearchEngine
    SearchEngine(std::vector<WebPage>* pages){
        this->pages = pages;
        std::srand(std::time(NULL));
        std::cout << (std::rand()/ (double) INT_MAX) << std::endl;
        websiteID = static_cast<int>((std::rand()/ (double) RAND_MAX) * pages->size());

        for(int i = 0; i < pages->size(); i++){
            webPageMap[pages->at(i).getWebsiteName()] = i;
        }
    }

    void startExecution(){
        while(true){
            std::string nextWebPage = pages->at(websiteID).loadWebPage();
            try{
                if(pages->at(websiteID).isWebsiteAvailable(nextWebPage)){
                    websiteID = webPageMap.at(nextWebPage);
                }else{
                    throw std::out_of_range("");
                }
            }catch(std::out_of_range){
                std::cout << "WebPage not found. Try again!" << std::endl;
            }
        }
    }
};

int main(){
    std::vector<double> poly1Vals{1, 0, -1};
    Polynomial poly1(poly1Vals);

    std::vector<double> poly2Vals{1, 0, -1};
    Polynomial poly2(poly2Vals);

    poly1.multiplyPolynomial(poly2);
    std::vector<double> zeros = poly1.getAllZeros();
    for(double d : zeros){
        std::cout << d << " is a zero of the fnc" << std::endl;
    }

    std::vector<double> polytVals{-1, -1, 1, 1};
    Polynomial polyt(polytVals);

    polyt.syntheticDivision(-1).printPolynomial();


    /*
    Rational rat = Rational{poly1, poly2};
    rat.printRational();
    rat.simplify();
    rat.printRational();
    */

}