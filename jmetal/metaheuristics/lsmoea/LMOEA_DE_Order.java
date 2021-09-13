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
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;
import jmetal.util.wrapper.XReal;

public class LMOEA_DE_Order extends Algorithm{
	
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
	private Operator crossoverOperator_;
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
	int gSize_;
	int numOfGroup_;
	int[] allVariables_;
	
	GroupingStrategy grouping;
	
	public LMOEA_DE_Order(Problem problem) {
		super(problem);
		numObj_ = problem.getNumberOfObjectives();
		numVar_ = problem.getNumberOfVariables();
	}
	
	public void initialization() throws JMException, ClassNotFoundException{
		generations_ = 0;
		populationSize_ = ((Integer) getInputParameter("populationSize")).intValue();
		//maxGenerations_ = ((Integer) getInputParameter("maxGenerations")).intValue();
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
		VectorGenerator vg1, vg2;
		vg1 = new TwoLevelWeightVectorGenerator(div1_, div2_,numObj_);
		lamada1_ = vg1.getVectors();
		
		if(numObj_ == 2){
			vg2 = new TwoLevelWeightVectorGenerator(99, 0,numObj_);
			maxGenerations_ = 5000;
		}else if(numObj_ == 3){
			vg2 = new TwoLevelWeightVectorGenerator(15, 0,numObj_);
			maxGenerations_ = 3334;
		}else{
			vg2 = new TwoLevelWeightVectorGenerator(5, 0,numObj_);
		}
		
		lamada2_ = vg2.getVectors();
		
		
		numOfGroup_ = 25;
		
		grouping = new GroupingStrategy(numOfGroup_, numVar_, numObj_);
		group_ = grouping.linearGrouping();
	
		allVariables_ = new int[numVar_];
		upper_ = new double[numVar_];
		lower_ = new double[numVar_];
		for(int i=0;i<numVar_;i++){
			allVariables_[i] = i;
			upper_[i] = problem_.getUpperLimit(i);
			lower_[i] = problem_.getLowerLimit(i);
		}
		
		//Create the initial population
		Solution newSolution;
		for (int i = 0; i < populationSize_; i++) {
			newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			population_.add(newSolution);
		} // for
	}

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		initialization();
		SolutionSet oldPopulation = new SolutionSet();
		while(generations_ < maxGenerations_){
			for(int i=0;i<populationSize_;i++){
				oldPopulation.add(population_.get(i));
			}
			reproduction();
			environmentalSelection();
			group_ = grouping.orderedGrouping(oldPopulation, population_);
			oldPopulation.clear();
			generations_++;
		}
		NondominatedRanking final_ranking = new NondominatedRanking(population_);
		return final_ranking.getSubfront(0);
	}
	
	public void reproduction() throws JMException{
		estimateIdealPoint(population_);
		normalizedPopulation(population_);
		SolutionSet[] solSet = new PopulationClassification(population_,numObj_,lamada2_).classification();
		int winnerSize = solSet[0].size();
		int loserSize = solSet[1].size();
		for(int i = 0; i < loserSize; i++){
			//Solution parent1 = solSet[1].get(i);
			Solution[] parent = new Solution[2];
			Solution child1 = new Solution(solSet[1].get(i));
			XReal xChild1 = new XReal(child1);
			double value;
			for(int j=0;j<numOfGroup_;j++){
				int rd1 = PseudoRandom.randInt(0, solSet[0].size()-1);
				int rd2 = PseudoRandom.randInt(0, solSet[0].size()-1);
				while(rd1 == rd2){
					rd2 = PseudoRandom.randInt(0, solSet[0].size()-1);
				}
				parent[0] = solSet[0].get(rd1);
				parent[1] = solSet[0].get(rd2);
				XReal xParent1 = new XReal(parent[0]);
				XReal xParent2 = new XReal(parent[1]);
				double val;
				for(int m=0;m<group_[j].length;m++){
					val = xChild1.getValue(group_[j][m]) //+ 0.5*(xParent1.getValue(group_[j][m]) - xChild1.getValue(group_[j][m]))
							+0.5*(xParent1.getValue(group_[j][m]) - xParent2.getValue(group_[j][m]));
					if (val < lower_[group_[j][m]])
						val = lower_[group_[j][m]];
					if (val > upper_[group_[j][m]])
						val = upper_[group_[j][m]];
					xChild1.setValue(group_[j][m], val);
				}
			}
			// Apply mutation
			mutationOperator_.execute(new Object[]{child1,group_});
			// Evaluation
			problem_.evaluate(child1);
			offspringPopulation_.add(child1);
		}
		
		for(int i = 0; i < winnerSize; i++){
			Solution child2 = new Solution(solSet[0].get(i));
			mutationOperator_.execute(new Object[]{child2,group_});
			// Evaluation
			problem_.evaluate(child2);
			offspringPopulation_.add(child2);
		}
	} 
	
	public void environmentalSelection(){
		// Create the solutionSet union of solutionSet and offSpring
		unionPopulation_ = ((SolutionSet) population_).union(offspringPopulation_);
		population_.clear();
		offspringPopulation_.clear();
		SolutionSet st = getStSolutionSet(unionPopulation_,populationSize_);
		estimateIdealPoint(st);
		normalizedPopulation(st);
		SDEDistanceAssign(st);
		population_ = new PopulationClassification(st,numObj_,lamada1_).classification()[0];
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
