#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <algorithm>

#include <shamir.h>

std::string lower_string(std::string &str)
{
    for (int i = 0; i < str.length(); i++)
    {
        str[i] = tolower(str[i]);
    }
    return str;
}

int main()
{
    int opt_number, opt_threshold;
    std::string secret;
    std::vector<std::string> shares = {};

    std::cout << "Enter the number of shares: ";
    std::cin >> opt_number;

    std::cout << "Enter the threshold: ";
    std::cin >> opt_threshold;

    std::cout << "Enter the secret:";
    std::cin >> secret;

    split(opt_threshold, opt_number, secret, shares);

    // Print the shares
    for (int i = 0; i < shares.size(); i++)
    {
        std::cout << "Share " << i + 1 << ": " << shares[i] << std::endl;
    }


    std::cout << "Trying to combine " << opt_threshold << " shares:" << std::endl;

    std::vector<bool> select(opt_number);
    std::fill(select.begin(), select.begin() + opt_threshold, true);

    do{
        std::vector<std::string> secret_shares;
        for(int i = 0; i < opt_number; ++i){
            if(select[i]){
                secret_shares.push_back(shares[i]);
                std::cout << i << ",";
            }
        }
        std::cout << ":";

        std::string combined_secret;
        combine(secret_shares, combined_secret);
        if(lower_string(secret) == lower_string(combined_secret))
        {
            std::cout << "The secret is successfully combined" << std::endl;
        }
        else
        {
            std::cout << "The secret is not successfully combined" << std::endl;
            exit(1);
        }
        std::cout << "The combined secret is: " << combined_secret << std::endl;
    }while(std::prev_permutation(select.begin(), select.end()));

    std::cout << std::endl;
}