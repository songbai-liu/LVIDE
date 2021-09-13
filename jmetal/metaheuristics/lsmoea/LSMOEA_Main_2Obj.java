//  LSMOEA_Main.java
//
// Author:
//     Songbai Liu <Songbai209@qq.com>
// Copyright (c) 2019 Songbai Liu
//
// This Program is free software: you can redistribute it and/or modify 
// it under the terms of the GNU Lesser General Public License as published
// by the free software foundation, either version 4 of license, or any 
// later version (at your option).

package jmetal.metaheuristics.lsmoea;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.nsgaII.NSGADE;
import jmetal.metaheuristics.nsgaII.NSGAII1;
import jmetal.metaheuristics.smpso.SMPSO;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.Fonseca;
import jmetal.problems.Kursawe;
import jmetal.problems.Schaffer;
import jmetal.problems.DTLZ.DTLZ1;
import jmetal.problems.DTLZ.DTLZ2;
import jmetal.problems.DTLZ.DTLZ3;
import jmetal.problems.DTLZ.DTLZ4;
import jmetal.problems.DTLZ.DTLZ5;
import jmetal.problems.DTLZ.DTLZ6;
import jmetal.problems.DTLZ.DTLZ7;
import jmetal.problems.IMF.IMF1;
import jmetal.problems.IMF.IMF10;
import jmetal.problems.IMF.IMF2;
import jmetal.problems.IMF.IMF3;
import jmetal.problems.IMF.IMF4;
import jmetal.problems.IMF.IMF5;
import jmetal.problems.IMF.IMF6;
import jmetal.problems.IMF.IMF7;
import jmetal.problems.IMF.IMF8;
import jmetal.problems.IMF.IMF9;
import jmetal.problems.LSMOP.LSMOP1;
import jmetal.problems.LSMOP.LSMOP2;
import jmetal.problems.LSMOP.LSMOP3;
import jmetal.problems.LSMOP.LSMOP4;
import jmetal.problems.LSMOP.LSMOP5;
import jmetal.problems.LSMOP.LSMOP6;
import jmetal.problems.LSMOP.LSMOP7;
import jmetal.problems.LSMOP.LSMOP8;
import jmetal.problems.LSMOP.LSMOP9;
import jmetal.problems.M2M.MOP1;
import jmetal.problems.M2M.MOP2;
import jmetal.problems.M2M.MOP3;
import jmetal.problems.M2M.MOP4;
import jmetal.problems.M2M.MOP5;
import jmetal.problems.M2M.MOP6;
import jmetal.problems.M2M.MOP7;
import jmetal.problems.MaF.MaF1;
import jmetal.problems.MaF.MaF13;
import jmetal.problems.MaF.MaF2;
import jmetal.problems.MaF.MaF3;
import jmetal.problems.MaF.MaF4;
import jmetal.problems.MaF.MaF5_Convex;
import jmetal.problems.MaF.MaF6;
import jmetal.problems.MaF.MaF7;
import jmetal.problems.MaF.MaF8;
import jmetal.problems.WFG.WFG1;
import jmetal.problems.WFG.WFG2;
import jmetal.problems.WFG.WFG3;
import jmetal.problems.WFG.WFG4;
import jmetal.problems.WFG.WFG5;
import jmetal.problems.WFG.WFG6;
import jmetal.problems.WFG.WFG7;
import jmetal.problems.WFG.WFG8;
import jmetal.problems.WFG.WFG9;
import jmetal.problems.ZDT.*;
import jmetal.problems.cec2009Competition.*;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfghvCalculator1;
import jmetal.util.Configuration;
import jmetal.util.JMException;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

/**
 * Class to configure and execute the NSGA-II algorithm.
 * 
 * Besides the classic NSGA-II, a steady-state version (ssNSGAII) is also
 * included (See: J.J. Durillo, A.J. Nebro, F. Luna and E. Alba "On the Effect
 * of the Steady-State Selection Scheme in Multi-Objective Genetic Algorithms"
 * 5th International Conference, EMO 2009, pp: 183-197. April 2009)
 */

public class LSMOEA_Main_2Obj {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	/**
	 * @param args
	 *            Command line arguments.
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException
	 *             Usage: three options -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main problemName -
	 *             jmetal.metaheuristics.nsgaII.NSGAII_main problemName
	 *             paretoFrontFile
	 */
	public static void printGD(String path,double[] GD){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path)     ;//java文件输出流，创建文件流
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;//OutputStreamWriter是字符流通向字节流的桥梁 
	      BufferedWriter bw      = new BufferedWriter(osw)        ;//缓冲区               
	      for (int i = 0; i < GD.length; i++) {  
	        bw.write(GD[i]+" ");//写到缓冲区
	        bw.newLine(); //换行       
	      }
	      
	      /* Close the file */
	      bw.close();
	    }catch (IOException e) {
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();
	    }       
	  } // printGD
	
	public static void printave(String path,double aveIGD,double varianceIGD,double aveHypervolume,double varianceHV){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path) ;
	      OutputStreamWriter osw = new OutputStreamWriter(fos) ;
	      BufferedWriter bw      = new BufferedWriter(osw)  ;            
	            
	     // for (int i = 0; i < IGD.length; i++) {  
	        
	        bw.write(aveIGD+" ");
	        bw.newLine(); 
	        bw.write(varianceIGD+" ");
	        bw.newLine();
	        bw.write(aveHypervolume+" ");
	        bw.newLine();  
	        bw.write(varianceHV+" ");
	        bw.newLine();
	        /* Close the file */
		      bw.close();
		    }catch (IOException e) {
		      Configuration.logger_.severe("Error acceding to the file");
		      e.printStackTrace();
		    }       
		  } // printave
	
	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException {
		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("NSGAII_main.log");
		logger_.addHandler(fileHandler_);
		int n = 100;
		for(int fun=1;fun<=18;fun++){
			int runtimes=1;
			double[] GDarray=new double[runtimes];
			double[] IGDarray=new double[runtimes];
			double[] spreadarray=new double[runtimes];
			double[] Hypervolume=new double[runtimes];
			long Execution_time=0;
			Problem problem=null; // The problem to solve
		    Algorithm algorithm; // The algorithm to use
		    Operator crossover1; // Crossover operator
		    Operator crossover2; // Crossover operator
		    Operator mutation1; // Mutation operator
		    Operator mutation2; // Mutation operator
		    Operator selection; // Selection operator

		    HashMap parameters; // Operator parameters

		    QualityIndicator indicators; // Object to get quality indicators

		    indicators = null;
		    if (args.length == 1) {
		    	Object[] params = { "Real" };
		    	problem = (new ProblemFactory()).getProblem(args[0], params);
		    } // if
		    else if (args.length == 2) {
		    	Object[] params = { "Real" };
		    	problem = (new ProblemFactory()).getProblem(args[0], params);
		    	indicators = new QualityIndicator(problem, args[1]);
		    } // if
		    else { // Default problem
		    	if(fun==1){
		    		problem = new LSMOP1("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP1_3D_10000.txt");
			    }
		    	if(fun==2){
		    		problem = new LSMOP2("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP2_3D_10000.txt");
			    }
		    	if(fun==3){
		    		problem = new LSMOP3("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP3_3D_10000.txt");
			    }
		    	if(fun==4){
		    		problem = new LSMOP4("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP4_3D_10000.txt");
			    }
		    	if(fun==5){
		    		problem = new LSMOP5("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP5_3D_10000.txt");
		    	}
		    	if(fun==6){
		    		problem = new LSMOP6("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP6_3D_10000.txt");
			    }
		    	if(fun==7){
		    		problem = new LSMOP7("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP7_3D_10000.txt");
			    }
		    	if(fun==8){
		    		problem = new LSMOP8("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP8_3D_10000.txt");
			    }
		    	if(fun==9){
		    		problem = new LSMOP9("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\LSMOP\\LSMOP9_3D_10000.txt");
			    }
		    	if(fun==10){
		    		problem = new IMF4("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF4_3D_10000.txt");
			    }
		    	if(fun==11){
		    		problem = new IMF8("Real",n,3);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF8_3D_10000.txt");
			    }
		    	if(fun==12){
		    		problem = new IMF1("Real",n,2);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF1_2D_5000.txt");
			    }
		    	if(fun==13){
		    		problem = new IMF2("Real",n,2);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF2_2D_5000.txt");
			    }
		    	if(fun==14){
		    		problem = new IMF3("Real",n,2);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF3_2D_5000.txt");
			    }
		    	if(fun==15){
		    		problem = new IMF5("Real",n,2);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF5_2D_5000.txt");
			    }
		    	if(fun==16){
		    		problem = new IMF6("Real",n,2);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF6_2D_5000.txt");
			    }
		    	if(fun==17){
		    		problem = new IMF7("Real",n,2);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF7_2D_5000.txt");
			    }
		    	if(fun==18){
		    		problem = new IMF9("Real",n,2);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF9_2D_5000.txt");
			    }
		    	if(fun==19){
		    		problem = new IMF10("Real",n,2);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\IMF\\IMF10_2D_5000.txt");
			    }
		    	if(fun==20){
		    		problem = new UF1("Real",n);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF1_500.txt");
			    }
		    	if(fun==21){
		    		problem = new UF2("Real",n);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF2_500.txt");
			    }
		    	if(fun==22){
		    		problem = new UF3("Real",n);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF3_500.txt");
			    }
		    	if(fun==23){
		    		problem = new UF4("Real",n);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF4_500.txt");
			    }
		    	if(fun==24){
		    		problem = new UF5("Real",n,2,0.1);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF5_21.txt");
			    }
		    	if(fun==25){
		    		problem = new UF6("Real",n,2,0.1);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF6_668.txt");
			    }
		    	if(fun==26){
		    		problem = new UF7("Real",n);
		    		indicators = new QualityIndicator(problem,"D:\\Matlab\\TruePF\\UF\\UF7_500.txt");
			    }
		    	if(fun==27){
		    		problem = new MOP1("Real",n);
		    		indicators = new QualityIndicator(problem, "D:\\Matlab\\TruePF\\MOP\\MOP1_2D_5000.txt") ;
			    }
		    	if(fun==28){
		    		problem = new MOP2("Real",n);
		    		indicators = new QualityIndicator(problem, "D:\\Matlab\\TruePF\\MOP\\MOP2_2D_5000.txt") ;
			    }
		    	if(fun==29){
		    		problem = new MOP3("Real",n);
		    		indicators = new QualityIndicator(problem, "D:\\Matlab\\TruePF\\MOP\\MOP3_2D_5000.txt") ;
			    }
		    	if(fun==30){
		    		problem = new MOP4("Real",n);
		    		indicators = new QualityIndicator(problem, "D:\\Matlab\\TruePF\\MOP\\MOP4_2D_5000.txt") ;
			    }
		    	if(fun==31){
		    		problem = new MOP5("Real",n);
		    		indicators = new QualityIndicator(problem, "D:\\Matlab\\TruePF\\MOP\\MOP5_2D_5000.txt") ;
			    }
			
		    } // else
		    for(int i=0;i<runtimes;i++){
		    	//algorithm = new LVIDE(problem);
		    	algorithm = new LVIDE_Analysis_Importance(problem);
		    	// Algorithm parameters
		    	if(fun <= 11){
		    		algorithm.setInputParameter("populationSize", 300);
		    		//algorithm.setInputParameter("maxEvaluations", Gmax);
		    		algorithm.setInputParameter("T", 30);
		    		algorithm.setInputParameter("div1", 23);
		    		algorithm.setInputParameter("div2", 0);
		    	}else{
		    		algorithm.setInputParameter("populationSize", 200);
		    		//algorithm.setInputParameter("maxEvaluations", Gmax);
		    		algorithm.setInputParameter("T", 20);
		    		algorithm.setInputParameter("div1", 199);
		    		algorithm.setInputParameter("div2", 0);
		    	}
		    	// Mutation and Crossover for Real codification
		    	// Crossover operator
		    	parameters = new HashMap();
		    	parameters.put("CR", 1.0);
		    	parameters.put("F", 0.5);
		    	crossover1 = CrossoverFactory.getCrossoverOperator("SegmentalDECrossover", parameters);
		    	
		    	parameters = new HashMap();
		    	parameters.put("probability", 1.0);
		    	parameters.put("distributionIndex", 30.0);
		    	crossover2 = CrossoverFactory.getCrossoverOperator("SegmentalSBXCrossover", parameters);
		    	
	    	 
		    	parameters = new HashMap();
		    	parameters.put("probability", 1.0 / problem.getNumberOfVariables());
		    	parameters.put("distributionIndex", 20.0);
		    	mutation1 = MutationFactory.getMutationOperator("GroupPolynomialMutation",parameters);//Group
		    	mutation2 = MutationFactory.getMutationOperator("PolynomialMutation",parameters);

		    	// Selection Operator
		    	parameters = null;
		    	selection = SelectionFactory.getSelectionOperator("RandomSelection",
		    			parameters);

		    	// Add the operators to the algorithm
		    	algorithm.addOperator("crossover1", crossover1);
		    	algorithm.addOperator("crossover2", crossover2);
		    	algorithm.addOperator("mutation1", mutation1);
		    	algorithm.addOperator("mutation2", mutation2);
		    	algorithm.addOperator("selection", selection);

		    	// Add the indicator object to the algorithm
		    	algorithm.setInputParameter("indicators", indicators);

		    	// Execute the Algorithm
		    	long initTime = System.currentTimeMillis();
		    	SolutionSet population = algorithm.execute();
		    	Execution_time+=(System.currentTimeMillis() - initTime);

		    	// Result messages
		    	//population.printObjectivesToFile("LVIDE_Rand_"+problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ "_" + problem.getNumberOfVariables() + "D_run"+(i+1)+".txt" );
		    	IGDarray[i]=indicators.getIGD1(population);
		    }
		    /* wfghvCalculator1 wfg = new wfghvCalculator1(population,fun);
		  		Hypervolume[i] = wfg.calculatewfghv();
			}
			printGD("NSGAII_10Obj_"+problem.getName()+"_HV.txt",Hypervolume);
			printGD("NSGAII_10Obj_"+problem.getName()+"_IGD.txt",IGDarray);*/
		    double sumIGD=0;
		    for(int i=0;i<runtimes;i++){
		    	sumIGD+=IGDarray[i];
		    }	  	  
		    System.out.println("IGD_"+ problem.getNumberOfObjectives()+"Obj_"+problem.getName()+ 
				  "_" + problem.getNumberOfVariables() + "D" + " = "+sumIGD/runtimes);
		} //main
	}
} // NSGAII_main
