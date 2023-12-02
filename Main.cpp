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
#include <string.h>





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

    static Matrix generateScalarMatrix(double scalar, int dimension){
        std::vector<std::vector<Term>> mat;

        for(int row = 0; row < dimension; row++){
            std::vector<Term> tempRow;
            for(int col = 0; col < dimension; col++){
                tempRow.push_back(Term(0, scalar));
            }
            mat.push_back(tempRow);
        }
        return Matrix(mat);
    }

    void scaleMatrix(double scalar){
        for(int row = 0; row < rows; row++){
            for(int col = 0; col < cols; col++){
                mat.at(row).at(col).setConstant(scalar * mat.at(row).at(col).getConstant());
            }
        }
    }
    

    void addMatrix(Matrix additionalMatrix){
        if(additionalMatrix.getCols() != cols) throw std::runtime_error("Matrices Cols are Different Sizes");

        if(additionalMatrix.getRows() != rows) throw std::runtime_error("Matrices Rows are Different Sizes");

        for(int row = 0; row < rows; row++){
            for(int col = 0; col < cols; col++){
                Term tempTerm(0, additionalMatrix.getTerm(row, col).getConstant() + this->getTerm(row, col).getConstant());
                setTerm(tempTerm, row, col);
            }
        }
    }

    static Matrix generateIdentity(double scalar, int dimension){
        std::vector<std::vector<Term>> mat;

        for(int row = 0; row < dimension; row++){
            std::vector<Term> tempRow;
            for(int col = 0; col < dimension; col++){
                if(col == row){
                    tempRow.push_back(Term(0, scalar));
                }else{
                    tempRow.push_back(Term(0, 0));
                }
            }
            mat.push_back(tempRow);
        }
        return Matrix(mat);
    }

    void putInRREF(){
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
                determinant *= pivotValue;
                //With the leading row determined we need to scale it to be 1
                for(int column = columnStart; column < cols; column++){
                    mat.at(leadingRow).at(column).setConstant(mat.at(leadingRow).at(column).getConstant()/pivotValue);
                }

                for(int row = 0; row < rows; row++){
                    if(row != leadingRow){
                        double leadingValue = mat.at(row).at(columnStart).getConstant();
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
    }

    std::vector<std::vector<double>> getKernelBasis(){
        std::vector<std::vector<double>> basisList;

        std::vector<int> pivotRowsIndexes;
        for(int col = 0; col < cols; col++){
            bool t = !isPivotColumn(col, &pivotRowsIndexes);
            if(t){
                std::vector<double> basisVector;
                for(int row = 0; row < rows; row++){
                    std::vector<Term> currRow = mat.at(row);
                    if(currRow.at(col).getConstant() == 0 && row == col){
                        basisVector.push_back(1);
                    }else{
                        basisVector.push_back(-currRow.at(col).getConstant());
                    }
                }
                for(double d : basisVector){
                    std::cout << d << std::endl;
                }
                basisList.push_back(basisVector);
            }
        }
        return basisList;
    }

    static void printBasisList(std::vector<std::vector<double>> kernelBasis){
        for(std::vector<double> vector : kernelBasis){
            for(double scalar : vector){
                std::cout << scalar << " ";
            }
            std::cout << std::endl;
        }
    }

    double calcDeterminant(){
        double determinant = 1;

        std::vector<std::vector<Term>> tempMat = mat;

        std::vector<int> pivotRows;
        if(rows != cols){
            throw std::runtime_error("The Matrix is not square");
        }

        //Find the pivot column
        for (int columnStart = 0; columnStart < cols; columnStart++){
            int leadingRow = -1;
            double pivotValue = 0;
            for(int row = 0; row < rows; row++){
                double value = tempMat.at(row).at(columnStart).getConstant();
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
                    tempMat.at(leadingRow).at(column).setConstant(tempMat.at(leadingRow).at(column).getConstant()/pivotValue);
                }

                for(int row = 0; row < rows; row++){
                    if(row != leadingRow){
                        double leadingValue = tempMat.at(row).at(columnStart).getConstant();
                        std::cout << leadingValue << "leadingValue" << std::endl;
                        for(int column = columnStart; column < cols; column++){
                            tempMat.at(row).at(column).setConstant(tempMat.at(row).at(column).getConstant() - leadingValue * tempMat.at(leadingRow).at(column).getConstant());
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

    int getRows(){
        return rows;
    }

    int getCols(){
        return cols;
    }

    void setTerm(Term term, int rowIndex, int colIndex){
        mat.at(rowIndex).at(colIndex) = term;
    }

    Term getTerm(int rowIndex, int colIndex){
        return mat.at(rowIndex).at(colIndex);
    }
    //Columns Start At 0
    bool isPivotColumn(int col, std::vector<int>* pivotRowsIndexes){
        int numOnes = 0;
        for(int rowIndex = 0; rowIndex < rows; rowIndex++){
            std::vector<Term> row = mat.at(rowIndex);
            if(row.at(col).getConstant() == 1 && !isVectorContains(*pivotRowsIndexes, rowIndex)){
                numOnes++;
                pivotRowsIndexes->push_back(rowIndex);
                if(numOnes > 1){
                    return false;
                }
            }else if(row.at(col).getConstant() != 0){
                return false;
            }
        }
        return true;
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
            zeros.push_back(zero);
        }

        return zeros;
    }

    Polynomial syntheticDivision(double divisor){
        std::vector<double> dividend;
        double value = polynomial.at(polynomial.size() - 1);
        dividend.push_back(value);
        for(int i = (polynomial.size() - 2); i >= 1; i--){
            value = value * divisor;
            value += polynomial.at(i);
            dividend.push_back(value);
        }
        for(double d : dividend){
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

    //Unique ID for the webpage
    int webID;

    public:
    //Constructor for the webpage
    WebPage(std::string websiteName, std::string keyword, int webID){
        this->websiteName = websiteName;
        this->keyword = keyword;
        this->webID =webID;
    }

    int getID(){
        return webID;
    }

    std::vector<double> getTransitionVector(std::vector<int> idOrder){
        std::vector<double> transitionVector;
        for(int id : idOrder){
            bool isFindLink = false;
            for(WebPage page : links){
                std::cout << "id " << id << " pageID " << page.getID() << std::endl;

                if(id == page.getID()){
                    transitionVector.push_back(((double) 1/links.size()));
                    isFindLink = true;
                }
            }
            if(!isFindLink){
                transitionVector.push_back(0);
            }
        }
        return transitionVector;
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

    std::string getKeyword(){
        return keyword;
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

    std::vector<WebPage> search(std::string search){
        std::vector<WebPage> orderedPages;
        std::vector<WebPage> searchedPages = getSearchedForPages(search);
        for(auto p : searchedPages){
            std::cout << p.getWebsiteName() <<std::endl;
        }
        std::vector<int> pageIDOrder = generatePageIDOrderVector(searchedPages);
        for(auto p : pageIDOrder){
            std::cout << p <<std::endl;
        }
        Matrix transitionMatrix = generateTransitionMatrix(pageIDOrder);
        transitionMatrix.addMatrix(Matrix::generateIdentity(-1, pageIDOrder.size()));
        transitionMatrix.putInRREF();
        std::vector<std::vector<double>> kernelBasis = transitionMatrix.getKernelBasis();
        if(kernelBasis.empty()) throw std::runtime_error("Kernel Empty");
        std::vector<double> equilibriumVector = kernelBasis[0];
        for(auto p : equilibriumVector){
            std::cout << p <<std::endl;
        }
        int eqSize = equilibriumVector.size();
        for(int m = 0; m < eqSize; m++){
            int maxIndex = 0;
            for(int i = 0; i < equilibriumVector.size(); i++){
                if(equilibriumVector.at(i) > equilibriumVector.at(maxIndex)){
                    maxIndex = i;
                }
            }
            orderedPages.push_back(searchedPages.at(maxIndex));
            searchedPages.erase(searchedPages.begin() + maxIndex);
            equilibriumVector.erase(equilibriumVector.begin() + maxIndex);
        }
        return orderedPages;
    }

    std::vector<WebPage> getSearchedForPages(std::string search){
        std::vector<WebPage> searchedPages;
        for(WebPage page : *pages){
            if(strcmp(search.c_str(), page.getKeyword().c_str()) == 0){
                searchedPages.push_back(page);
            }
        }
        return searchedPages;
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

    std::vector<int> generatePageIDOrderVector(std::vector<WebPage> pages){
        //ID order for the webpages
        std::vector<int> pageIDOrder;
        for(WebPage page : pages){
            pageIDOrder.push_back(page.getID());
        }
        return pageIDOrder;
    }

    Matrix generateTransitionMatrix(std::vector<int> pageIDOrder){
        std::vector<std::vector<Term>> transitionMatrix;
        for(int i = 0; i < pages->size(); i++){
            std::vector<Term> row;
            for(int j = 0; j < pages->size(); j++){
                row.push_back(Term(0, 0));
            }
            transitionMatrix.push_back(row);
        }
        for(int col = 0; col < pages->size(); col++){
            WebPage page = pages->at(col);
            std::vector<double> transitionVector = page.getTransitionVector(pageIDOrder);
            for(int row = 0; row < pageIDOrder.size(); row++){
                transitionMatrix.at(row).at(col) = Term(0, transitionVector.at(row));
            }
        }
            
        return Matrix(transitionMatrix);
    }
};

int main(){
    WebPage google("Google", "Search Engine", 0);
    WebPage bing("Bing", "Search Engine", 1);
    WebPage duckDuckGo("DuckDuck", "Search Engine", 2);

    //Set Up Google's Links
    std::vector<WebPage> googleLinks{bing, duckDuckGo};
    google.setLinks(googleLinks);

    //Set Up Bing's Links
    std::vector<WebPage> bingLinks{google, duckDuckGo};
    bing.setLinks(bingLinks);

    //Set Up Duck Duck Go's Links
    std::vector<WebPage> duckDuckGoLinks{bing};
    duckDuckGo.setLinks(duckDuckGoLinks);

    //Set up the Engine Over the Entire Search Space
    std::vector<WebPage> allLinks{google, bing, duckDuckGo};
    SearchEngine engine(&allLinks);

    std::cout << "Before gen" << std::endl;
    auto pages = engine.search("Search Engine");
    for(WebPage page : pages){
        std::cout << page.getWebsiteName() << std::endl;
    }

    Matrix square = Matrix::generateScalarMatrix(1,3);
    square.printMatrix();
    square.scaleMatrix(0.5);
    square.printMatrix();
    //engine.startExecution();
}