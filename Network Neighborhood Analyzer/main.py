"""
Author: P.A.R.R. Wijerathne
Date: 3/26/2024
Project 7: Network neighborhood analyzer - Implementation(Main Program)

"""

from NetworkAnalyzer import *
import os

print("Welcome to Network Analyzer")
print("---------------------------")

print("Enter the required details to continue")
while True:
    pff = input("ENTER THE NAME OF THE PROTEIN FUNCTION FILE IN THE DIRECTORY: ")
    if not os.path.isfile(pff):
        print("---GIVEN FILE NAME DOES NOT EXIST IN THE DIRECTORY---")
    else:
        break

while True:
    ppi = input("ENTER THE NAME OF THE TSV FILE CONTAINING PPI NETWORK INFORMATION: ")
    if not os.path.isfile(ppi):
        print("---GIVEN FILE NAME DOES NOT EXIST IN THE DIRECTORY---")
    else:
        break

print("Choose an option from below")
print("----------------------------")
print("1. THE NUMBER OF PROTEINS AND INTERACTIONS IN THE GIVEN PPI NETWORK")
print("2. FIND THE DEGREE OF A SPECIFIC PROTEIN")
print("3. GET THE DEGREE DISTRIBUTION OF A PPI NETWORK")
print("4. MAJORITY VOTING - GET THE COUNT OF PROTEINS ANNOTATED TO THE GIVEN FUNCTION IN THE IMMEDIATE NEIGHBORHOOD")
print("5. ANALYZE THE NEIGHBORHOOD OF A SPECIFIC PROTEIN")
print("6. DISTRIBUTION OF THE COUNT OF PROTEINS ANNOTATED TO THE GIVEN FUNCTION IN THE IMMEDIATE NEIGHBORHOOD")
print("7. HISHIGHAKI METHOD - PREDICT CANDIDATE GENES FOR A SINGLE FUNCTION USING THE HISHIGAKI ALGORITHM")
print("8. HISTOGRAM OF HISHIGAKI SCORE DISTRIBUTION")
print("9. BOXPLOT OF HISHIGAKI SCORE DISTRIBUTION")
print("10. EXIT")

analyzer = NetworkAnalyzer(ppi)

while True:
    print("-------------------------------------------------------------------")
    num = int(input("Enter the number relevant for the functionality: "))

    if num==1:
        print(analyzer.count_proteins_interactions())
    elif num ==2:
        protein = input("Enter the protein: ")
        print(analyzer.calculate_degree(protein))
    elif num==3:
        analyzer.plot_degree_distribution()
    elif num==4:
        print(analyzer.neighborhood_pcounter_order(pff))
    elif num==5:
        protein = input("Enter the protein: ")
        analyzer.protein_neighbors(pff, protein)
    elif num==6:
        print(analyzer.annotated_protein_distribution(pff, True))
    elif num==7:
        print(analyzer.hishigaki_analyze(pff))
    elif num ==8:
        analyzer.plot_hishigaki_score_distribution(pff)
    elif num ==9:
        analyzer.hishigaki_box_plot(pff)
    elif num ==10:
        break
    else:
        print("Input a valid number")

