
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

// setup functions to help us along 

void create_cell_types( void );
void setup_tissue( void ); 

// set up the BioFVM microenvironment 
void setup_microenvironment( void ); 

// custom pathology coloring function 

std::vector<std::string> my_coloring_function( Cell* );

// custom functions can go here 

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt );
void custom_function( Cell* pCell, Phenotype& phenotype , double dt );

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt ); 

void stem_cell_rule();
void divide_custom_data();
int cancer_cell_rule( Cell* );
