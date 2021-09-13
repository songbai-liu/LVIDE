package jmetal.metaheuristics.lsmoea;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.moead.Utils;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;
import jmetal.util.wrapper.XReal;

public class LMOEA1 extends Algorithm{
	
	private int populationSize_;//population size
	private int numObj_; //number of objectives
	private int numVar_; //number of variables
	/**
	 * Stores the population
	 */
	private SolutionSet population_;
	private SolutionSet offspringPopulation_;
	private SolutionSet unionPopulation_;
	
	private int evaluations_;
	private int maxEvaluations_;
	
	/**
	 * Operators
	 */
	private Operator crossoverOperator_;
	private Operator mutationOperator_;
	private Operator selectionOperator_;
	
	double[] upper_;
	double[] lower_;
	
	double[] zideal_;
	double[][] lamada_;
	int div1_,div2_;
	
	private SolutionSet detection_Set_;
	double[][] improvement_;
	int[][] group_;
	int gSize_;
	int numOfGroup_;
	
	public LMOEA1(Problem problem) {
		super(problem);
		numObj_ = problem.getNumberOfObjectives();
		numVar_ = problem.getNumberOfVariables();
	}
	
	public void initialization() throws JMException, ClassNotFoundException{
		evaluations_ = 0;
		populationSize_ = ((Integer) getInputParameter("populationSize")).intValue();
		maxEvaluations_ = ((Integer) getInputParameter("maxEvaluations")).intValue();
		div1_ = ((Integer) getInputParameter("div1")).intValue();
		div2_ = ((Integer) getInputParameter("div2")).intValue();
		population_ = new SolutionSet(populationSize_);
		offspringPopulation_ = new SolutionSet(populationSize_);
		unionPopulation_ = new SolutionSet(2*populationSize_);
		
		//Read the operators
		mutationOperator_ = operators_.get("mutation");
		crossoverOperator_ = operators_.get("crossover");	
		selectionOperator_ = operators_.get("selection");
		zideal_ = new double[numObj_];
		/* generate two-layer weight vectors */
		VectorGenerator vg;
		vg = new TwoLevelWeightVectorGenerator(div1_, div2_,numObj_);
		lamada_ = vg.getVectors();
		
		upper_ = new double[numVar_];
		lower_ = new double[numVar_];
		
		detection_Set_ = new SolutionSet();
		improvement_ = new double[numVar_][numObj_];
		numOfGroup_ = 10;
		group_ = new int[numOfGroup_][];
		gSize_ = numVar_/numOfGroup_;
		for(int g=0;g<numOfGroup_;g++){
			group_[g] = new int[gSize_];
		}
		
		variablesDetection();
		adaptiveGrouping();
		//detection_Set_.printObjectivesToFile("Detection_"+numObj_+"Obj_"+problem_.getName()+ "_" + numVar_ + "D_run" + ".txt" );
		estimateIdealPoint(detection_Set_);
		normalizedPopulation(detection_Set_);
		//Create the initial population
		for (int i = 0; i < populationSize_; i++) {
			population_ = new PopulationClassification(detection_Set_,numObj_,lamada_).classification()[0];
		} // for
	}

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		initialization();
		int initEvaluations = evaluations_;
		int T_ = 0;
		int[] inGroup = group_[0];
		while(evaluations_ < maxEvaluations_){
			reproduction(inGroup);
			environmentalSelection();
			if(T_ < numOfGroup_-1){
				if((evaluations_ - initEvaluations)%70000 == 0){
					T_++;
					int[] groupA = inGroup;
					int[] groupB = group_[T_];
					inGroup = new int[groupA.length + groupB.length];
					System.arraycopy(groupA, 0, inGroup, 0, groupA.length);
					System.arraycopy(groupB, 0, inGroup, groupA.length, groupB.length);
				}
			}
		}
		NondominatedRanking final_ranking = new NondominatedRanking(population_);
		return final_ranking.getSubfront(0);
	}
	
	public void reproduction(int[] group) throws JMException{
		Solution[] parents = new Solution[2];
		for (int i = 0; i < (populationSize_)/2; i++) {
			// obtain parents
			parents = (Solution[]) selectionOperator_.execute(population_);
			//parents[1] = (Solution) selectionOperator_.execute(population_);
			Solution[] offSpring = (Solution[]) crossoverOperator_.execute(new Object[]{parents,group});
			//Solution[] offSpring = (Solution[]) crossoverOperator_.execute(parents);
			mutationOperator_.execute(new Object[]{offSpring[0],group_});
			mutationOperator_.execute(new Object[]{offSpring[1],group_});
			//mutationOperator_.execute(offSpring[0]);
			//mutationOperator_.execute(offSpring[1]);
			problem_.evaluate(offSpring[0]);
			problem_.evaluateConstraints(offSpring[0]);
			problem_.evaluate(offSpring[1]);
			problem_.evaluateConstraints(offSpring[1]);
			offspringPopulation_.add(offSpring[0]);
			offspringPopulation_.add(offSpring[1]);
			evaluations_ += 2;	
		} // for
	}
	
	public void environmentalSelection(){
		// Create the solutionSet union of solutionSet and offSpring
		unionPopulation_ = ((SolutionSet) population_).union(offspringPopulation_);
		population_.clear();
		offspringPopulation_.clear();
		SolutionSet st = getStSolutionSet(unionPopulation_,populationSize_);
		estimateIdealPoint(st);
		normalizedPopulation(st);
		population_ = new PopulationClassification(st,numObj_,lamada_).classification()[0];
	}
	
	/*
	 * Estimate the Ideal Point 
	 */
	public void estimateIdealPoint(SolutionSet solutionSet){
		for(int i=0; i<numObj_;i++){
			zideal_[i] = 1.0e+30;
			for(int j=0; j<solutionSet.size();j++){
				if(solutionSet.get(j).getObjective(i) < zideal_[i]){
					zideal_[i] = solutionSet.get(j).getObjective(i);
				}//if
			}//for
		}//for
	}
	
	public void normalizedPopulation(SolutionSet solutionSet){
		double minNomal = Double.MAX_VALUE;
		int minIndex = 0;
		for(int i=0; i<solutionSet.size();i++){
			Solution sol = solutionSet.get(i);
			double sum = 0.0;
			double normal = 0.0;
			for(int j=0; j<problem_.getNumberOfObjectives();j++){
				double value = sol.getObjective(j) -  zideal_[j];
				sol.setNormalizedObjective(j, value);
				sum += value;
				normal += value*value;
			}
			normal = Math.sqrt(normal);
			if(normal < minNomal) {
				minNomal = normal;
				minIndex = i;
			}
			sol.setDistanceToIdealPoint(normal);
			sol.setSumValue(sum);
		}
	}
	
	public SolutionSet getStSolutionSet(SolutionSet ss,int size) {
		Ranking ranking = new NondominatedRanking(ss);
		int remain = size;
		int index = 0;
		SolutionSet front = null;
		SolutionSet mgPopulation = new SolutionSet();
		front = ranking.getSubfront(index);
		while ((remain > 0) && (remain >= front.size())) {

			for (int k = 0; k < front.size(); k++) {
				mgPopulation.add(front.get(k));
			} // for
			// Decrement remain
			remain = remain - front.size();
			// Obtain the next front
			index++;
			if (remain > 0) {
				front = ranking.getSubfront(index);
			} // if
		}
		if (remain > 0) { // front contains individuals to insert
			for (int k = 0; k < front.size(); k++) {
				mgPopulation.add(front.get(k));
			}
		}
		return mgPopulation;
	}
	
	public void variablesDetection() throws ClassNotFoundException, JMException{
		//int numVar_ = problem_.getNumberOfVariables();
		//int numObj_ = problem_.getNumberOfObjectives();
		int N_ = numObj_; //number of solutions for detection
		int k = 20; //number of segments for detection on each variable
		double[][] orders = new double[numVar_][k+1];//record the segments
		double[] detalD = new double[numVar_];
		for(int d=0;d<numVar_;d++){
			detalD[d] = (problem_.getUpperLimit(d) - problem_.getLowerLimit(d))/4;
			double interval = (problem_.getUpperLimit(d) - problem_.getLowerLimit(d))/k;//the length of each segment
			for(int j=0;j<k;j++){
				orders[d][j] = problem_.getLowerLimit(d) + j*interval;
			}
			orders[d][k] = problem_.getUpperLimit(d);
		}//for
		
		double[][][] improvement = new double[N_][numVar_][numObj_];
		
		Solution[] sol = new Solution[N_];//solutions for detection
		XReal[] individual = new XReal[N_];
		
		int[][] oIndex = new int[N_][numVar_];//record the original order index
		
		for(int n=0;n<N_;n++){
			sol[n] = new Solution(problem_);
			individual[n] = new XReal(sol[n]);
			for(int i=0; i<numVar_;i++){
				int rd = PseudoRandom.randInt(0, k);
				oIndex[n][i] = rd;
				individual[n].setValue(i, orders[i][rd]);
			}//for
			problem_.evaluate(sol[n]);
			//evaluations_++;
			detection_Set_.add(new Solution(sol[n]));
		}
		
		double[][] subFit = new double[N_][numObj_];
		double[][] bestSubFit = new double[N_][numObj_];
		
		
		for(int n=0;n<N_;n++){
			for(int m=0;m<numObj_;m++){
				subFit[n][m] = sol[n].getObjective(m);
				bestSubFit[n][m] = sol[n].getObjective(m);
			}
		}
		
		int[][][] bestLevel = new int[N_][numObj_][numVar_];
		int[][][] worstLevel = new int[N_][numObj_][numVar_];
		double[][][] minFitness = new double[N_][numObj_][numVar_];
		double[][][] maxFitness = new double[N_][numObj_][numVar_];
		
		//double[][][] improvement = new double[N_][numVar_][numObj_];
		
		double fit = 0;
		double[] oObjective = new double[numObj_];
		int minIndex;
		/*
		  get the best position of the order and the average contribution
		  for each variable on each subproblem
		*/
		for(int n=0;n<N_;n++){
			for(int d=0;d<numVar_;d++){
				for(int m=0;m<numObj_;m++){
					bestLevel[n][m][d] = oIndex[n][d];
					worstLevel[n][m][d] = oIndex[n][d];
					minFitness[n][m][d] = sol[n].getObjective(m);
					maxFitness[n][m][d] = sol[n].getObjective(m);
					oObjective[m] = sol[n].getObjective(m);
				}//m
				for(int j=0; j<=k; j++){
					if(j != oIndex[n][d]){
						individual[n].setValue(d, orders[d][j]);
						problem_.evaluate(sol[n]);
						detection_Set_.add(new Solution(sol[n]));
						evaluations_++;
						for(int m=0;m<numObj_;m++){
							fit = sol[n].getObjective(m);
							if(fit < minFitness[n][m][d]){
								minFitness[n][m][d] = fit;
								bestLevel[n][m][d] = j;
							}
							if(fit > maxFitness[n][m][d]){
								maxFitness[n][m][d] = fit;
								worstLevel[n][m][d] = j;
							}
						}
					}
				}//j
				
				for(int m=0;m<numObj_;m++){
					improvement[n][d][m] = oObjective[m] - minFitness[n][m][d];
					//improvement[n][d][m] = maxFitness[n][m][d] - minFitness[n][m][d];
				}
				individual[n].setValue(d, orders[d][bestLevel[n][n][d]]);
				problem_.evaluate(sol[n]);
			}//d	
		}//n
		
		for(int d=0;d<numVar_;d++){
			for(int m=0;m<numObj_;m++){
				improvement_[d][m] = improvement[m][d][m];
			}
		}	
	}//variable detection

	public void adaptiveGrouping(){
		//int numVar_ = problem_.getNumberOfVariables();
		//int numObj_ = problem_.getNumberOfObjectives();
		// dominateMe[i] contains the number of solutions dominating i
		int[] dominateMe = new int[numVar_];
		// iDominate[k] contains the list of solutions dominated by k
		List<Integer>[] iDominate = new List[numVar_];
		// front[i] contains the list of individuals belonging to the front i
		List<Integer>[] front = new List[numVar_ + 1];
		// flagDominate is an auxiliar encodings.variable
		int flagDominate;
		
		// Initialize the fronts
		for (int i = 0; i < front.length; i++)
			front[i] = new LinkedList<Integer>();
		for (int d = 0; d < numVar_; d++) {
			// Initialize the list of individuals that i dominate and the number
			// of individuals that dominate me
			iDominate[d] = new LinkedList<Integer>();
			dominateMe[d] = 0;
		}
		
		for (int p = 0; p < (numVar_ - 1); p++) {
			for (int q = p + 1; q < numVar_; q++) {
				flagDominate = dominance_Compare(improvement_[p],improvement_[q]);
				if (flagDominate == -1) {
					iDominate[p].add(q);
					dominateMe[q]++;
				} else if (flagDominate == 1) {
					iDominate[q].add(p);
					dominateMe[p]++;
				}
			}
		}
		
		for (int p = 0; p < numVar_; p++) {
			if (dominateMe[p] == 0) {
				front[0].add(p);
			}
		}
		
		// Obtain the rest of fronts
		int i = 0;
		Iterator<Integer> it1, it2; // Iterators
		while (front[i].size() != 0) {
			i++;
			it1 = front[i - 1].iterator();
			while (it1.hasNext()) {
				it2 = iDominate[it1.next()].iterator();
				while (it2.hasNext()) {
					int index = it2.next();
					dominateMe[index]--;
					if (dominateMe[index] == 0) {
						front[i].add(index);
					}
				}
			}
		}
		int[][] group = new int[i][];
		for(int j=0;j<i;j++){
			group[j] = new int[front[j].size()];
			for(int f=0;f<front[j].size();f++){
				group[j][f] = front[j].get(f);
			}
		}
		
		//Contribution-based grouping
		int L = group.length;
		int[] tt = new int[numVar_];
		int u = 0;
		for(int l=0;l<L;l++){
			for(int t=0;t<group[l].length;t++){
				tt[u] = group[l][t];
				u++;
			}
		}
		
		int t = 0;
		for(int g=0;g<numOfGroup_;g++){
			for(int m=0;m<gSize_;m++){
				group_[g][m] = tt[t];
				t++;
			}
		}

	}
	
	public int dominance_Compare(double[] value1, double[] value2){
		//int numObj_ = problem_.getNumberOfObjectives();
		if (value1==null)
		    return 1;
		else if (value2 == null)
		    return -1;

		int dominate1 ; // dominate1 indicates if some objective of solution1 
		                    // dominates the same objective in solution2. dominate2
		int dominate2 ; // is the complementary of dominate1.

		dominate1 = 0 ; 
		dominate2 = 0 ;
		    
		int flag; //stores the result of the comparison
		    
		// Equal number of violated constraints. Applying a dominance Test then
		double fit1, fit2;
		for (int i = 0; i < numObj_; i++) {
		   fit1 = value1[i];
		   fit2 = value2[i];
		   if (fit1 > fit2) {
		      flag = -1;
		   } else if (fit1 < fit2) {
		      flag = 1;
		   } else {
		      flag = 0;
		   }
		      
		   if (flag == -1) {
		      dominate1 = 1;
		   }
		      
		   if (flag == 1) {
		      dominate2 = 1;           
		   }
		}
		            
		if (dominate1 == dominate2) {            
		   return 0; //No one dominate the other
		}
		if (dominate1 == 1) {
		   return -1; // solution1 dominate
		}
		return 1;    // solution2 dominate  
	}//dominance_Compare

}
