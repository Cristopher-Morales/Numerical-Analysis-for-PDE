import runpy

def problems_assignment_1(current_directory, problem_to_run):
	if problem_to_run==2:
		runpy.run_path(current_directory+"/assignment1/Problem2/problem2numericalanalysis_report2.py")
		print("successfully executed")
	elif problem_to_run==3:
		runpy.run_path(current_directory+"/assignment1/Problem3/problem3numericalanalysis_report2.py")
		print("successfully executed")
	elif problem_to_run==4:
		runpy.run_path(current_directory+"/assignment1/Problem4/problem4numericalanalysis_report2.py")
		print("successfully executed")
	else:
		print("That is not a valid problem to be run.Please try again")
		print("Unsuccessfully executed")
def problems_assignment_2(current_directory, problem_to_run):
	if problem_to_run==1:
		runpy.run_path(current_directory+"/assignment2/Problem1/problem1numericalanalysis_report2.py")
		print("successfully executed")
	elif problem_to_run==2:
		runpy.run_path(current_directory+"/assignment2/Problem2/problem2numericalanalysis_report2.py")
		print("successfully executed")
	elif problem_to_run==3:
		runpy.run_path(current_directory+"/assignment2/Problem3/problem3numericalanalysis_report2.py")
		print("successfully executed")
	elif problem_to_run==4:
		runpy.run_path(current_directory+"/assignment2/Problem4/problem4numericalanalysis_report2.py")
		print("successfully executed")
	elif problem_to_run==5:
		runpy.run_path(current_directory+"/assignment2/Problem5/problem5numericalanalysis_report2.py")
		print("successfully executed")	
	else: 
		print("That is not a valid problem to be run.Please try again")
		print("Unsuccessfully executed")
