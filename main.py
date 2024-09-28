import runpy
import os
import shutil
import assignments as asg

def main():
    current_directory = os.getcwd()
    assignment_to_run=int(input("enter number of the assignment to run (for example, if you want to run assigment 1, please insert 1:"))
    if assignment_to_run==1:
        problem_to_run=int(input("enter number of problem to run (for example, if you want to run problem 1, please insert 1:"))
        asg.problems_assignment_1(current_directory, problem_to_run)
    elif assignment_to_run==2:
        problem_to_run=int(input("enter number of problem to run (for example, if you want to run problem 1, please insert 1:"))
        asg.problems_assignment_2(current_directory, problem_to_run)	
    else:
        print("That is not a valid assignment to be run.Please try again")
        print("Unsuccessfully executed")
    shutil.rmtree(current_directory+"/__pycache__")

if __name__ == "__main__":
    main()
