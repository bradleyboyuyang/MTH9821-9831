/**
 * test program for our ProductServices
 */

#include <iostream>
#include "products.hpp"
#include "productservice.hpp"

using namespace std;

// Color codes
#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN    "\033[36m"
#define WHITE   "\033[37m"

int main()
{
  std::cout << BLUE << "=================================Testing Bond Class=================================" << RESET << std::endl;

  // Create the 10Y treasury note
  date maturityDate(2025, Nov, 16);
  string cusip = "912828M56";
  Bond treasuryBond(cusip, CUSIP, "T", 2.25, maturityDate);

  // Create the 2Y treasury note
  date maturityDate2(2017, Nov, 5);
  string cusip2 = "912828TW0";
  Bond treasuryBond2(cusip2, CUSIP, "T", 0.75, maturityDate2);

  // Create a BondProductService
  BondProductService *bondProductService = new BondProductService();

  // Add the 10Y note to the BondProductService and retrieve it from the service
  bondProductService->Add(treasuryBond);
  Bond bond = bondProductService->GetData(cusip);
  cout << "CUSIP: " << bond.GetProductId() << " ==> " << bond << endl;

  // Add the 2Y note to the BondProductService and retrieve it from the service
  bondProductService->Add(treasuryBond2);
  bond = bondProductService->GetData(cusip2);
  cout << "CUSIP: " << bond.GetProductId() << " ==> " << bond << endl;

  std::cout << BLUE << "=============================Testing Interest Rate Swap=============================" << RESET << std::endl;
  // Create the Spot 10Y Outright Swap
  date effectiveDate(2015, Nov, 16);
  date terminationDate(2025, Nov, 16);
  string outright10Y = "Spot-Outright-10Y";
  IRSwap outright10YSwap(outright10Y, THIRTY_THREE_SIXTY, THIRTY_THREE_SIXTY, SEMI_ANNUAL, LIBOR, TENOR_3M, effectiveDate, terminationDate, USD, 10, SPOT, OUTRIGHT);

  // Create the IMM 2Y Outright Swap
  date effectiveDate2(2015, Dec, 20);
  date terminationDate2(2017, Dec, 20);
  string imm2Y = "IMM-Outright-2Y";
  IRSwap imm2YSwap(imm2Y, THIRTY_THREE_SIXTY, THIRTY_THREE_SIXTY, SEMI_ANNUAL, LIBOR, TENOR_3M, effectiveDate2, terminationDate2, USD, 2, IMM, OUTRIGHT);

  // Create a IRSwapProductService
  IRSwapProductService *swapProductService = new IRSwapProductService();

  // Add the Spot 10Y Outright Swap to the IRSwapProductService and retrieve it from the service
  swapProductService->Add(outright10YSwap);
  IRSwap swap = swapProductService->GetData(outright10Y);
  cout << "Swap: " << swap.GetProductId() << " ==> " << swap << endl;

  // Add the IMM 2Y Outright Swap to the IRSwapProductService and retrieve it from the service
  swapProductService->Add(imm2YSwap);
  swap = swapProductService->GetData(imm2Y);
  cout << "Swap: " << swap.GetProductId() << " ==> " << swap << endl;

  // Question 2: test three futures: FX, Eurodollar, Bond
  std::cout << BLUE << "===================================Testing Futures===================================" << RESET << std::endl;
  // test FX Future
  date futureContractDate(2016, Dec, 16);
  string futureContractID = "FX12345";
  Future fxFuture(futureContractID, FX, NONE_TYPE, CME, DEC, 0.01, "FX_Base", 100, futureContractDate);
  cout << "Future: " << fxFuture.GetProductId() << " ==> " << fxFuture << endl;

  // Test Eurodollar Future
  float liborRate = 2.5;
  date euroFutureDate(2019, Oct, 11);
  string euroFutureID = "Euro12345";
  EurodollarFuture eurodollarFuture(euroFutureID, INTEREST_RATE, EURODOLLAR, CME, DEC, 0.005, "Eurodollar_Base", 100, euroFutureDate, liborRate);
  cout << "Eurodollar Future: " << eurodollarFuture.GetProductId() << " ==> " << eurodollarFuture << endl;

  // Test Bond Future
  BondFutureType bondFutureType = TEN_YR;
  string bondFutureID = "Bond12345";
  date bondFutureDate(2023, Jan, 1);
  BondFuture bondFuture(bondFutureID, INTEREST_RATE, BOND_FUTURE, CME, MAR, 0.01, "Bond_Base", 100, bondFutureDate, bondFutureType);
  cout << "Bond Future: " << bondFuture.GetProductId() << " ==> " << bondFuture << endl;

  // Test FutureProductService and retrive three futures
  std::cout << BLUE << "==================================Retriving Futures==================================" << RESET << std::endl;
  FutureProductService *futureProductService = new FutureProductService();
  futureProductService->Add(fxFuture);
  futureProductService->Add(eurodollarFuture);
  futureProductService->Add(bondFuture);
  Future future = futureProductService->GetData(futureContractID);
  cout << "Future: " << future.GetProductId() << " ==> " << future << endl;
  future = futureProductService->GetData(euroFutureID);
  cout << "Eurodollar Future: " << future.GetProductId() << " ==> " << future << endl;
  future = futureProductService->GetData(bondFutureID);
  cout << "Bond Future: " << future.GetProductId() << " ==> " << future << endl;

  // Question 3: Test GetBonds and GetSwaps Functions
  std::cout << BLUE << "===========================Test GetBonds and GetSwaps Functions==========================" << RESET << std::endl;
  // Test GetBonds Function
  string ticker = "T";
  vector<Bond> bonds = bondProductService->GetBonds(ticker);
  cout << CYAN << "Bonds with ticker " << ticker << " are: " << RESET << endl;
  for (auto it = bonds.begin(); it != bonds.end(); ++it)
  {
    cout << *it << endl;
  }

  // Test GetSwaps with the specified fixed leg day count convention
  DayCountConvention fixedLegDayCountConvention = THIRTY_THREE_SIXTY;
  vector<IRSwap> swaps = swapProductService->GetSwaps(fixedLegDayCountConvention);
  cout << CYAN << "Swaps with fixed leg day count convention 30/360"  << " are: " << RESET << endl;
  for (auto it = swaps.begin(); it != swaps.end(); ++it)
  {
    cout << *it << endl;
  }

  // Test GetSwaps with specified fixed leg payment frequency
  PaymentFrequency fixedLegPaymentFrequency = SEMI_ANNUAL;
  swaps = swapProductService->GetSwaps(fixedLegPaymentFrequency);
  cout << CYAN << "Swaps with fixed leg payment frequency Semi-Annual" << " are: " << RESET << endl;
  for (auto it = swaps.begin(); it != swaps.end(); ++it)
  {
    cout << *it << endl;
  }

  // Test GetSwaps with the specified floating index
  FloatingIndex floatingIndex = LIBOR;
  swaps = swapProductService->GetSwaps(floatingIndex);
  cout << CYAN << "Swaps with floating index LIBOR" << " are: " << RESET << endl;

  for (auto it = swaps.begin(); it != swaps.end(); ++it)
  {
    cout << *it << endl;
  }

  // Test GetSwaps with a term in years greater than the specified value
  int termInYears = 5;
  swaps = swapProductService->GetSwapsGreaterThan(termInYears);
  cout << CYAN << "Swaps with a term in years greater than " << termInYears << " are: " << RESET << endl;
  for (auto it = swaps.begin(); it != swaps.end(); ++it)
  {
    cout << *it << endl;
  }

  // Test GetSwaps with the specified swap type
  SwapType swapType = IMM;
  swaps = swapProductService->GetSwaps(swapType);
  cout << CYAN << "Swaps with swap type IMM" << " are: " << RESET << endl;
  for (auto it = swaps.begin(); it != swaps.end(); ++it)
  {
    cout << *it << endl;
  }

  // Test GetSwaps with the specified swap leg type
  SwapLegType swapLegType = OUTRIGHT;
  swaps = swapProductService->GetSwaps(swapLegType);
  cout << CYAN << "Swaps with swap leg type OUTRIGHT" << " are: " << RESET << endl;
  for (auto it = swaps.begin(); it != swaps.end(); ++it)
  {
    cout << *it << endl;
  }

  return 0;

}
