package jmetal.metaheuristics.lsmoea;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.JMException;
import jmetal.util.Permutation;
import jmetal.util.PseudoRandom;
import jmetal.util.ranking.NondominatedRanking;
import jmetal.util.ranking.Ranking;
import jmetal.util.ranking.ThetaRanking;
import jmetal.util.vector.TwoLevelWeightVectorGenerator;
import jmetal.util.vector.VectorGenerator;
import jmetal.util.wrapper.XReal;

public class LSMOEA extends Algorithm{
	
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
	
	double[] zideal_;
	double[][] lamada1_;
	double[][] lamada2_;
	int div1_,div2_;
	int div3_,div4_;
	
	double[] upper_;
	double[] lower_;
	double[] interval_;
	double[] average_;
	int K_;

	public LSMOEA(Problem problem) {
		super(problem);
		numObj_ = problem.getNumberOfObjectives();
		numVar_ = problem.getNumberOfVariables();
	}

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		initialization();
		while(evaluations_ < maxEvaluations_){
			reproduction();
			environmentalSelection();
		}
		NondominatedRanking final_ranking = new NondominatedRanking(population_);
		return final_ranking.getSubfront(0);
	}
	
	public void initialization() throws JMException, ClassNotFoundException{
		evaluations_ = 0;
		populationSize_ = ((Integer) getInputParameter("populationSize")).intValue();
		maxEvaluations_ = ((Integer) getInputParameter("maxEvaluations")).intValue();
		div1_ = ((Integer) getInputParameter("div1")).intValue();
		div2_ = ((Integer) getInputParameter("div2")).intValue();
		div3_ = ((Integer) getInputParameter("div3")).intValue();
		div4_ = ((Integer) getInputParameter("div4")).intValue();
		population_ = new SolutionSet(populationSize_);
		offspringPopulation_ = new SolutionSet(populationSize_);
		unionPopulation_ = new SolutionSet(2*populationSize_);
		
		//Read the operators
		mutationOperator_ = operators_.get("mutation");
		crossoverOperator_ = operators_.get("crossover");
		
		//Create the initial population
		Solution newSolution;
		for (int i = 0; i < populationSize_; i++) {
			newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			population_.add(newSolution);
		} // for
		
		zideal_ = new double[problem_.getNumberOfObjectives()];
		/* generate two-layer weight vectors */
		VectorGenerator vg;
		vg = new TwoLevelWeightVectorGenerator(div1_, div2_,numObj_);
		lamada1_ = vg.getVectors();
		
		vg = new TwoLevelWeightVectorGenerator(div3_, div4_,numObj_);
		lamada2_ = vg.getVectors();
		upper_ = new double[numVar_];
		lower_ = new double[numVar_];
		interval_ = new double[numVar_];
		average_ = new double[numVar_];
		K_ = 50;
		for(int d=0;d<numVar_;d++){
			upper_[d] = problem_.getUpperLimit(d);
			lower_[d] = problem_.getLowerLimit(d);
			interval_[d] = (upper_[d] - lower_[d])/K_;
		}
	}
	
	public void reproduction() throws JMException{
		Ranking ranking = new NondominatedRanking(population_);
		estimateIdealPoint(population_);
		normalizedPopulation(population_);
		SolutionSet[] solSet = new PopulationClassification(population_,numObj_,lamada2_).classification();
		if(ranking.getSubfront(0).size() > population_.size()/2){
			updateAveragePopint(ranking.getSubfront(0));
			for (int i = 0; i < populationSize_/2; i++) {
				Solution parent1 = solSet[0].get(i);
				Solution parent2 = solSet[1].get(i);
				Solution child1 = new Solution(parent1);
				XReal xChild1 = new XReal(child1);
				double value;
				double rate = 1.0;
				for(int var=0;var<numVar_;var++){
					if(xChild1.getValue(var) > average_[var]){
						value = xChild1.getValue(var) + rate*PseudoRandom.randDouble(2*interval_[var],3*interval_[var]);
					}else{
						value = xChild1.getValue(var) + rate*PseudoRandom.randDouble(-3*interval_[var],-2*interval_[var]);
					}
					if (value < xChild1.getLowerBound(var))
						value = xChild1.getLowerBound(var);
					if (value > xChild1.getUpperBound(var))
						value = xChild1.getUpperBound(var);
					xChild1.setValue(var, value);
				}
				
				Solution child2 = new Solution(parent2);
				XReal xChild2 = new XReal(child2);
				rate = 1.0;
				double val;
				for(int var=0;var<numVar_;var++){
					if(xChild2.getValue(var) > average_[var]){
						val = xChild2.getValue(var) + rate*PseudoRandom.
								randDouble(average_[var]-xChild2.getValue(var)-interval_[var],interval_[var]);
					}else{
						val = xChild2.getValue(var) + rate*PseudoRandom.
								randDouble(-interval_[var], average_[var]-xChild2.getValue(var)+interval_[var]);
					}
					if (val < xChild2.getLowerBound(var))
						val = xChild2.getLowerBound(var);
					if (val > xChild2.getUpperBound(var))
						val = xChild2.getUpperBound(var);
					xChild2.setValue(var, val);
				}

				// Apply mutation
				mutationOperator_.execute(child1);
				
				mutationOperator_.execute(child2);

				// Evaluation
				problem_.evaluate(child1);
				problem_.evaluate(child2);
				offspringPopulation_.add(child1);
				offspringPopulation_.add(child2);
				evaluations_ += 2;
			}
		}else{
			for (int i = 0; i < populationSize_/2; i++) {
				
				Solution parent1 = solSet[0].get(i);
				Solution parent2 = solSet[1].get(i);
				Solution parent;
				
				int rd1 = PseudoRandom.randInt(0, solSet[0].size()-1);
				int rd2 = PseudoRandom.randInt(0, solSet[1].size()-1);
				/*while(rd1 == rd2){
					rd2 = PseudoRandom.randInt(0, solSet[0].size()-1);
				}*/
				Solution child1 = new Solution(parent1);
				XReal xChild1 = new XReal(child1);
				XReal xParent1 = new XReal(solSet[0].get(rd1));
				XReal xParent2 = new XReal(solSet[1].get(rd2));
				double value;
				for(int var=0;var<numVar_;var++){
					value = xChild1.getValue(var) + 0.5*(xParent1.getValue(var) - xParent2.getValue(var));
					if (value < xChild1.getLowerBound(var))
						value = xChild1.getLowerBound(var);
					if (value > xChild1.getUpperBound(var))
						value = xChild1.getUpperBound(var);
					xChild1.setValue(var, value);
				}
				
				Solution child2 = new Solution(parent2);
				XReal xChild2 = new XReal(child2);
				int rd = PseudoRandom.randInt(0, solSet[0].size()-1);;
				
				parent = solSet[0].get(rd);
				XReal xParent = new XReal(parent);
				for(int var=0;var<numVar_;var++){
					double val = xChild2.getValue(var) + 0.5*(xParent.getValue(var) - xChild2.getValue(var));
					if (val < xChild2.getLowerBound(var))
						val = xChild2.getLowerBound(var);
					if (val > xChild2.getUpperBound(var))
						val = xChild2.getUpperBound(var);
					xChild2.setValue(var, val);
				}

				// Apply mutation
				mutationOperator_.execute(child1);
				
				mutationOperator_.execute(child2);

				// Evaluation
				problem_.evaluate(child1);
				problem_.evaluate(child2);
				offspringPopulation_.add(child1);
				offspringPopulation_.add(child2);
				evaluations_ += 2;
			}
		}

	} 
	
	public void updateAveragePopint(SolutionSet solSet) throws JMException{
		int solSize = solSet.size();
		XReal[] sols = new XReal[solSize];
		for(int i=0;i<solSize;i++ ){
			sols[i] = new XReal(solSet.get(i));
		}
		for(int var=0;var<numVar_;var++){
			average_[var] = 0;
			for(int i=0;i<solSize;i++){
				average_[var] += sols[i].getValue(var);
			}
			average_[var] = average_[var]/solSize;
		}
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
	
	public void environmentalSelection(){
		// Create the solutionSet union of solutionSet and offSpring
		unionPopulation_ = ((SolutionSet) population_).union(offspringPopulation_);
		population_.clear();
		offspringPopulation_.clear();
		SolutionSet st = getStSolutionSet(unionPopulation_,populationSize_);
		estimateIdealPoint(st);
		normalizedPopulation(st);
		population_ = new PopulationClassification(st,numObj_,lamada1_).classification()[0];
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
}
