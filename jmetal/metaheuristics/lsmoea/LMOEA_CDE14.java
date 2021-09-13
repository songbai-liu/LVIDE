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
import jmetal.util.DistanceShifted;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.comparators.CrowdingDistanceComparator;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;
import jmetal.util.wrapper.XReal;

public class LMOEA_CDE14 extends Algorithm{
	
	private int populationSize_;//population size
	private int numObj_; //number of objectives
	private int numVar_; //number of variables
	/**
	 * Stores the population
	 */
	private SolutionSet population_;
	private SolutionSet offspringPopulation_;
	private SolutionSet unionPopulation_;
	
	private int generations_;
	private int maxGenerations_;
	
	/**
	 * Operators
	 */
	private Operator crossoverOperator1_;
	private Operator crossoverOperator2_;
	private Operator mutationOperator_;
	private Operator selectionOperator_;
	
	double[] upper_;
	double[] lower_;
	
	double[] zideal_;
	double[][] lamada1_;
	double[][] lamada2_;
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
	int[][] broadenGroup_;
	int numOfGroup_;
	int evaluations_;
	
	public LMOEA_CDE14(Problem problem) {
		super(problem);
		numObj_ = problem.getNumberOfObjectives();
		numVar_ = problem.getNumberOfVariables();
	}
	
	public void initialization() throws JMException, ClassNotFoundException{
		generations_ = 0;
		populationSize_ = ((Integer) getInputParameter("populationSize")).intValue();
		//maxGenerations_ = ((Integer) getInputParameter("maxEvaluations")).intValue();
		T_ = ((Integer) this.getInputParameter("T")).intValue();
		div1_ = ((Integer) getInputParameter("div1")).intValue();
		div2_ = ((Integer) getInputParameter("div2")).intValue();
		population_ = new SolutionSet(populationSize_);
		offspringPopulation_ = new SolutionSet(populationSize_);
		unionPopulation_ = new SolutionSet(2*populationSize_);
		
		neighborhood_ = new int[populationSize_][T_];
		
		//Read the operators
		mutationOperator_ = operators_.get("mutation1");
		crossoverOperator1_ = operators_.get("crossover1");
		crossoverOperator2_ = operators_.get("crossover2");	
		selectionOperator_ = operators_.get("selection");
		zideal_ = new double[numObj_];
		/* generate two-layer weight vectors */
		VectorGenerator vg1,vg2;
		vg1 = new TwoLevelWeightVectorGenerator(div1_, div2_,numObj_);
		lamada1_ = vg1.getVectors();
		
		if(numObj_ == 2){//N=200
			vg2 = new TwoLevelWeightVectorGenerator(99, 0,numObj_);
			maxGenerations_ = 5000;
		}else if(numObj_ == 3){//N = 300
			vg2 = new TwoLevelWeightVectorGenerator(14, 6,numObj_);
			maxGenerations_ = 3334;
		}else{
			vg2 = new TwoLevelWeightVectorGenerator(5, 0,numObj_);
			maxGenerations_ = 10;
		}
		lamada2_ = vg2.getVectors();
		
		
		upper_ = new double[numVar_];
		lower_ = new double[numVar_];
		
		detection_Set_ = new SolutionSet();
		improvement_ = new double[numVar_][numObj_];
		
		numOfGroup_ = numVar_/25;

		for(int i=0;i<numVar_;i++){
			upper_[i] = problem_.getUpperLimit(i);
			lower_[i] = problem_.getLowerLimit(i);
		}
		evaluations_ = 0;
		//variablesDetection();
		maxGenerations_ = maxGenerations_ - evaluations_/populationSize_;
		GroupingStrategy groupingStrategy = new GroupingStrategy(numOfGroup_,numVar_,numObj_);
		//group_ = groupingStrategy.contributionGrouping(improvement_);
		//memorySet_ = new SolutionSet[numOfGroup_];
		group_ = groupingStrategy.linearGrouping();
		//group_ = groupingStrategy.randomGrouping();
		broadenGroup_ = groupingStrategy.broadeningGrouping(group_);
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
		SDEDistanceAssign(detection_Set_);
		population_ = new PopulationClassification(detection_Set_,numObj_,lamada1_).classification()[0];
	}

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		initialization();
		while(generations_ <= maxGenerations_){
			reproduction_GSBX(broadenGroup_[numOfGroup_-1]);
			environmentalSelection();
			generations_++;
		}
		NondominatedRanking final_ranking = new NondominatedRanking(population_);
		return final_ranking.getSubfront(0);
	}
	
	public void reproduction_GDE(int[] group) throws JMException{
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
			matingSelection(p, n, 3, type);
			parents[0] = population_.get(p.get(0));
			parents[1] = population_.get(p.get(1));
			//parents[0] = population_.get(n);
			parents[2] = population_.get(p.get(2));
		
			//STEP2. Reproduction: Crossover + Mutation
			Solution offSpring = (Solution) crossoverOperator1_.execute(new Object[]{population_.get(n), parents,group});
			mutationOperator_.execute(new Object[]{offSpring,group_});
			//mutationOperator_.execute(offSpring);
			
			//STEP3. Function Evaluation
			problem_.evaluate(offSpring);
			problem_.evaluateConstraints(offSpring);
			offspringPopulation_.add(offSpring);
			//evaluations_ ++;	
		} // for
	}
	
	public void reproduction_GSBX(int[] group) throws JMException{
		estimateIdealPoint(population_);
		normalizedPopulation(population_);
		getNeighborhood();

		Solution[] parents = new Solution[2];
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
		
			//STEP2. Reproduction: Crossover + Mutation
			Solution[] offSpring = (Solution[]) crossoverOperator2_.execute(new Object[]{parents,group});
			mutationOperator_.execute(new Object[]{offSpring[0],group_});
			//mutationOperator1_.execute(offSpring[0]);
			
			//STEP3. Function Evaluation
			problem_.evaluate(offSpring[0]);
			problem_.evaluateConstraints(offSpring[0]);
			offspringPopulation_.add(offSpring[0]);
			//evaluations_ ++;	
		} // for
	}
	
	public void environmentalSelection(){
		// Create the solutionSet union of solutionSet and offSpring
		unionPopulation_ = ((SolutionSet) population_).union(offspringPopulation_);
		population_.clear();
		offspringPopulation_.clear();
		estimateIdealPoint(unionPopulation_);
		normalizedPopulation(unionPopulation_);
		SDEDistanceAssign(unionPopulation_);
		population_ = new PopulationClassification(unionPopulation_,numObj_,lamada1_).classification()[0];
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
	
	public void variablesDetection() throws ClassNotFoundException, JMException{
		int N_ = numObj_; //number of solutions for detection
		int k = 10; //number of segments for detection on each variable
		if(numVar_ == 100){
			k = 25;
		}else if(numVar_ == 300){
			k = 25;
		}else if(numVar_ == 500){
			k = 25;
		}else if(numVar_ == 800){
			k = 25;
		}else if(numVar_ == 1000){
			k = 20;
		}else if(numVar_ == 2000){
			k = 18;
		}else{
			k = 10;
		}
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
			evaluations_++;
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
						//detection_Set_.add(new Solution(sol));
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
				if(bestLevel == 0){
					individual.setValue(d, orders[d][bestLevel] + (0.01*detalD[d]/k));
				}else if(bestLevel == k){
					individual.setValue(d, orders[d][bestLevel] - (0.01*detalD[d]/k));
				}else{
					individual.setValue(d, orders[d][bestLevel] + (0.001*detalD[d]/k));
				}
				//individual.setValue(d, orders[d][bestLevel]);
				problem_.evaluate(sol);
				detection_Set_.add(new Solution(sol));
			}//d	
		}//n
	}//variable detection
	
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
  		while(list.size() < size){
  			if(type == 1){
  				r = PseudoRandom.randInt(0, ss - 1);
  				p = neighborhood_[cid][r];
  			}else{
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
  	
  	private void SDEDistanceAssign(SolutionSet winnerPopulation) {
		DistanceShifted distSDE = new DistanceShifted();
		for(int i =0; i < winnerPopulation.size(); i ++){
			Solution sol = winnerPopulation.get(i); 
			double disSDE = Double.MAX_VALUE;
			for(int j=0; j < i ; j++){ 
				Solution so2 = winnerPopulation.get(j);
				double tempDis = distSDE.distanceBetweenNormalizedObjectivesShifted(sol, so2);
				if(tempDis  < disSDE){
					disSDE = tempDis;
				} 
			}   
			sol.setCrowdingDistance(disSDE);
		}		
	} 

}
