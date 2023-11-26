#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <unordered_map>


struct WebPage{
    private:
    //Search Term
    std::string keyword;

    //Name of the Website
    std::string websiteName;

    //Vector of the WebPages this WebPage Links to
    std::vector<WebPage> links;

    public:
    //Constructor for the webpage
    WebPage(std::string websiteName, std::string keyword){
        this->websiteName = websiteName;
        this->keyword = keyword;
    }

    void setLinks(std::vector<WebPage> links){
        this->links = links;
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
};

class SearchEngine{
    private:
    //vector of the webpages inside this search engine
    std::vector<WebPage> pages;

    //Current WebPage
    int websiteID;

    //Map of WebPage names to their index
    std::unordered_map<std::string, int> webPageMap;

    public:
    //Constructor for the SearchEngine
    SearchEngine(std::vector<WebPage> pages){
        this->pages = pages;
        std::srand(std::time(nullptr));
        websiteID = static_cast<int>((std::rand()/INT_MAX) * pages.size());

        for(int i = 0; i < pages.size(); i++){
            webPageMap[pages.at(i).getWebsiteName()] = i;
        }
    }

    void startExecution(){
        while(true){
            std::string nextWebPage = pages.at(websiteID).loadWebPage();
            try{
                websiteID = webPageMap.at(nextWebPage);
            }catch(std::out_of_range){
                std::cout << "WebPage not found. Try again!" << std::endl;
            }
        }
    }
};

int main(){
    WebPage google("Google", "Google");
    WebPage bing("bing", "bing");
    std::vector<WebPage> googleLinks{bing, bing};
    google.setLinks(googleLinks);
    SearchEngine engine(googleLinks);
    engine.startExecution();
}