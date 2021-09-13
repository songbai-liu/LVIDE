package jmetal.metaheuristics.lsmoea;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;

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

public class LMOEA_DE1 extends Algorithm{
	
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
	
	/**
	 * T: neighbour size
	 */
	int T_;
	/**
	 * Neighborhood
	 */
	int[][] neighborhood_;
	
	private SolutionSet detection_Set_;
	double[][] improvement_;
	int[][] group_;
	int gSize_;
	int numOfGroup_;
	int[] allVariables_;
	
	public LMOEA_DE1(Problem problem) {
		super(problem);
		numObj_ = problem.getNumberOfObjectives();
		numVar_ = problem.getNumberOfVariables();
	}
	
	public void initialization() throws JMException, ClassNotFoundException{
		evaluations_ = 0;
		populationSize_ = ((Integer) getInputParameter("populationSize")).intValue();
		maxEvaluations_ = ((Integer) getInputParameter("maxEvaluations")).intValue();
		T_ = ((Integer) this.getInputParameter("T")).intValue();
		div1_ = ((Integer) getInputParameter("div1")).intValue();
		div2_ = ((Integer) getInputParameter("div2")).intValue();
		population_ = new SolutionSet(populationSize_);
		offspringPopulation_ = new SolutionSet(populationSize_);
		unionPopulation_ = new SolutionSet(2*populationSize_);
		
		neighborhood_ = new int[populationSize_][T_];
		
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
		numOfGroup_ = 25;
		group_ = new int[numOfGroup_][];
		gSize_ = numVar_/numOfGroup_;
		for(int g=0;g<numOfGroup_;g++){
			group_[g] = new int[gSize_];
		}
		allVariables_ = new int[numVar_];
		for(int i=0;i<numVar_;i++){
			allVariables_[i] = i;
			upper_[i] = problem_.getUpperLimit(i);
			lower_[i] = problem_.getLowerLimit(i);
		}
		
		variablesDetection();
		adaptiveGrouping();
		//detection_Set_.printObjectivesToFile("Detection_"+numObj_+"Obj_"+problem_.getName()+ "_" + numVar_ + "D_run" + ".txt" );
		
		//Create the initial population
		Solution newSolution;
		for (int i = 0; i < populationSize_; i++) {
			newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			detection_Set_.add(newSolution);
		} // for
		estimateIdealPoint(detection_Set_);
		normalizedPopulation(detection_Set_);
		population_ = new PopulationClassification(detection_Set_,numObj_,lamada_).classification()[0];
		
	}

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		initialization();
		int initEvaluations = evaluations_;
		int T_ = 0;
		int[] inGroup = group_[0];
		while(evaluations_ < maxEvaluations_){
			reproduction(inGroup);
			if(T_ < numOfGroup_-1){
				if((evaluations_ - initEvaluations)%30000 == 0){
					T_++;
					int[] groupA = inGroup;
					int[] groupB = group_[T_];
					inGroup = new int[groupA.length + groupB.length];
					System.arraycopy(groupA, 0, inGroup, 0, groupA.length);
					System.arraycopy(groupB, 0, inGroup, groupA.length, groupB.length);
				}
			}
			//reproduction(allVariables_);
			environmentalSelection();
		}
		NondominatedRanking final_ranking = new NondominatedRanking(population_);
		return final_ranking.getSubfront(0);
	}
	
	public void reproduction(int[] group) throws JMException{
		estimateIdealPoint(population_);
		normalizedPopulation(population_);
		getNeighborhood();
		Solution[] parents = new Solution[3];
		int[] permutation = new int[populationSize_];
		Utils.randomPermutation(permutation, populationSize_);
		for (int i = 0; i < populationSize_; i++) {
			int n = permutation[i];// int n = i ; // or int n = i;
			// obtain parents
			int type;
			double rnd = PseudoRandom.randDouble();

			// STEP1. Mating selection based on probability
			if (rnd < 0.8) // if (rnd < realb)
			{
				type = 1; // neighborhood
			} else {
				type = 2; // whole population
			}
			Vector<Integer> p = new Vector<Integer>();
			matingSelection(p, n, 2, type);
			parents[0] = population_.get(p.get(0));
			parents[1] = population_.get(p.get(1));
			parents[2] = population_.get(n);
		
			//STEP2. Reproduction: Crossover + Mutation
			Solution offSpring = (Solution) crossoverOperator_.execute(new Object[]{population_.get(n), parents,group});
			mutationOperator_.execute(new Object[]{offSpring,group_});
			//mutationOperator_.execute(offSpring);
			
			//STEP3. Function Evaluation
			problem_.evaluate(offSpring);
			problem_.evaluateConstraints(offSpring);
			offspringPopulation_.add(offSpring);
			evaluations_ ++;	
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
		int N_ = numObj_; //number of solutions for detection
		int k = 20; //number of segments for detection on each variable
		double[][] orders = new double[numVar_][k+1];//record the segments
		double[] detalD = new double[numVar_];
		for(int d=0;d<numVar_;d++){
			detalD[d] = (upper_[d] - lower_[d])/4;
			double interval = (upper_[d] - lower_[d])/k;//the length of each segment
			for(int j=0;j<k;j++){
				orders[d][j] = lower_[d] + j*interval;
			}
			orders[d][k] = upper_[d];
		}//for
		
		//double[][] improvement = new double[numVar_][numObj_];
		
		int[][] oIndex = new int[numVar_][numObj_];//record the original order index
		/*
		  get the best position of the order and the average contribution
		  for each variable on each subproblem
		*/
		for(int n=0;n<numObj_;n++){
			Solution sol = new Solution(problem_);
			XReal individual = new XReal(sol);
			for(int i=0; i<numVar_;i++){
				int rd = PseudoRandom.randInt(0, k);
				oIndex[i][n] = rd;
				individual.setValue(i, orders[i][rd]);
			}//for
			problem_.evaluate(sol);
			//evaluations_++;
			detection_Set_.add(new Solution(sol));

			for(int d=0;d<numVar_;d++){
				double oObjective = sol.getObjective(n);
				int bestLevel = oIndex[d][n];
				int worstLevel = oIndex[d][n];
				double minFitness = sol.getObjective(n);
				double maxFitness = sol.getObjective(n);
	
				for(int j=0; j<=k; j++){
					if(j != oIndex[d][n]){
						individual.setValue(d, orders[d][j]);
						problem_.evaluate(sol);
						detection_Set_.add(new Solution(sol));
						evaluations_++;
						
						double fit = sol.getObjective(n);
						if(fit < minFitness){
							minFitness = fit;
							bestLevel = j;
						}
						if(fit > maxFitness){
							maxFitness = fit;
							worstLevel = j;
						}
					}
				}//perform k times perturbations for each solution on each variable
				
				//improvement_[d][n] = oObjective - minFitness;
				improvement_[d][n] = maxFitness - minFitness;
				
				//individual.setValue(d, orders[d][worstLevel]);
				if(bestLevel == 0){
					individual.setValue(d, orders[d][bestLevel] + (0.5*detalD[d]/k));
				}else if(bestLevel == k){
					individual.setValue(d, orders[d][bestLevel] - (0.5*detalD[d]/k));
				}else{
					individual.setValue(d, orders[d][bestLevel]);
				}
				problem_.evaluate(sol);
				//sol = new Solution(sol);
				//individual = new XReal(sol);
			}//d	
		}//n
	}//variable detection

	public void adaptiveGrouping(){
		
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
		
		//Contribution-based grouping
		int L = front.length;
		int[] tt = new int[numVar_];
		int u = 0;
		for(int l=0;l<L;l++){
			for(int t=0;t<front[l].size();t++){
				tt[u] = front[l].get(t);
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
	
	/**
     * 
     */
  	public void getNeighborhood() {
  		double[] x = new double[populationSize_];
  		int[] idx = new int[populationSize_];

  		for (int i = 0; i < populationSize_; i++) {
  			// calculate the distances based on Solutions
  			for (int j = 0; j < populationSize_; j++) {
  				x[j] = computeAngle(population_.get(i), population_.get(j));
  				idx[j] = j;
  			} // for

  			// find 'niche' nearest neighboring solutions
  			Utils.minFastSort(x, idx, populationSize_, T_);
  			System.arraycopy(idx, 0, neighborhood_[i], 0, T_);
  		} // for
  	} // initNeighborhood
  	
  	 /*
     * Compute the angle value between Solution1 and Solution2
     */
	public double computeAngle(Solution s1, Solution s2){
		double angle = 0.0;
		double distanceToidealPoint1 = s1.getDistanceToIdealPoint();
		double distanceToidealPoint2 = s2.getDistanceToIdealPoint();
		double innerProduc = 0.0; 
		for(int i=0; i<problem_.getNumberOfObjectives(); i++){
			innerProduc += s1.getNormalizedObjective(i) * s2.getNormalizedObjective(i);
		}
		double value = innerProduc/(distanceToidealPoint1*distanceToidealPoint2);
		if(value > 1.0){
			value = 1.0;
		}
		angle = Math.acos(Math.abs(value));
		//System.out.println(Math.abs(innerProduc/(distanceToidealPoint1*distanceToidealPoint2)));
		return angle;
	}//computeAngle
	
	/**
     * 
     */
  	public void matingSelection(Vector<Integer> list, int cid, int size,
  			int type) {
  		// list : the set of the indexes of selected mating parents
  		// cid : the id of current subproblem
  		// size : the number of selected mating parents
  		// type : 1 - neighborhood; otherwise - whole population
  		int ss;
  		int r;
  		int p;
  		ss = neighborhood_[cid].length;
  		while (list.size() < size) {
  			if (type == 1) {
  				r = PseudoRandom.randInt(0, ss - 1);
  				p = neighborhood_[cid][r];
  			} else {
  				p = PseudoRandom.randInt(0, populationSize_ - 1);
  			}
  			boolean flag = true;
  			for (int i = 0; i < list.size(); i++) {
  				if (list.get(i) == p) // p is in the list
  				{
  					flag = false;
  					break;
  				}
  			}

  			// if (flag) list.push_back(p);
  			if (flag) {
  				list.addElement(p);
  			}
  		}
  	} // matingSelection

}
