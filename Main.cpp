#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <unordered_map>
#include <chrono>
#include <cmath>
#include <algorithm>

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
            currentPolynomial = currentPolynomial.syntheticDivision(zero);
            zeros.push_back(zero);
        }

        return zeros;
    }

    Polynomial syntheticDivision(double divisor){
        std::vector<double> dividend;
        int value = polynomial.at(polynomial.size() - 1);
        dividend.push_back(value);
        for(int i = (polynomial.size() - 2); i >= 1; i--){
            value *= divisor;
            value += polynomial.at(i);
            dividend.push_back(value);
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
    /*
    WebPage google("Google", "Search Engine");
    WebPage bing("Bing", "Search Engine");
    WebPage duckDuckGo("DuckDuck", "Search Engine");

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
    engine.startExecution();
    */
    std::vector<double> polynomialVals{6, -5, 1};
    Polynomial polynomial(polynomialVals);
    std::vector<double> zeros = {polynomial.getAllZeros()};
    for(double zero : zeros){
        std::cout << zero << std::endl;
    }
}