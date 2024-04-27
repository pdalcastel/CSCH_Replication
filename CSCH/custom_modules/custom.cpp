
#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	*/
        Custom_Cell_Data ccd;
        ccd.add_variable( "child_flag", "dimensionless", 0 );
        //ccd.add_variable( "division_capacity", "dimensionless", 10.0  );
        //std::cout << ccd << std::endl;

	/*
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 
///////////////////////////////////////


///////////////////////////////////////
void stem_cell_rule(void)
{ 
        // we want to take ccd.cancer_odds times monster_bias and apply to cell division to see if daughter cell is cancerous

	// simple check for life each turn
        static int stem_type = get_cell_definition( "stem cell" ).type;

        for( int i=0; i < (*all_cells).size() ;i++ )
        {
                Cell* pCell;
                pCell = (*all_cells)[i];

		if (pCell->phenotype.flagged_for_division == true) {
                	// let us know the cell is alive and flagged for division
			std::cout << "Stem cell " << pCell->ID << " of type " << pCell->type << " is alive and flagged for division" << std::endl;
        	}
        }

        //std::cout << ccd << std::endl;
        return;
}
///////////////////////////////////////


///////////////////////////////////////
int cancer_cell_rule(Cell* tested_cell)
{
	//display_cell_definitions( std::cout );
	static int gen_limit = 10;     //TODO use the XML config setting and not just be hard-coded
	static int sc_type = get_cell_definition( "stem cell" ).type; // returns index 0 for stem cell and 1 for cancer cell
	std::string sc_name = get_cell_definition( "stem cell" ).name; // returns index 0 for stem cell and 1 for cancer cell
	static int cc_type = get_cell_definition( "cancer" ).type; // returns index 0 for stem cell and 1 for cancer cell
	std::string cc_name = get_cell_definition( "cancer" ).name; // returns index 0 for stem cell and 1 for cancer cell
	//std::cout << "Stem cell: " << sc_name << " (" << sc_type << ")" << std::endl;
	//std::cout << "Cancer cell: " << cc_name << " (" << cc_type << ")" << std::endl;
	//std::cout << "tested_cell: " << tested_cell->type_name << " (" << tested_cell->type << ")" << std::endl;

	std::cout << "===== cancer_cell_rule is applied here to cell " << tested_cell->ID << " (type " << tested_cell->type << ") with generation count: " << tested_cell->custom_data[2] << " =====" << std::endl;

	//impose a limit on the number of times a cancer cell may divide before dying
	static int doom = 0;
	if(tested_cell->type == cc_type && tested_cell->custom_data[2] == gen_limit){
		//std::cout << "Cancer cell has reached the end of its allowed generation_count" << std::endl;
		doom = 0;
	} else if(tested_cell->type == cc_type) {
		//std::cout << "Cancer cell still viable" << std::endl;
		//tested_cell->custom_data[2] += 1;
		//std::cout << " - Cancer cell generation value incremented to: " << tested_cell->custom_data[2] << std::endl;
		doom = 0;
	} else {
		//tested_cell->custom_data[2] += 1;
		//std::cout << " - Stem cell generation value incremented to: " << tested_cell->custom_data[2] << std::endl;
		doom = 0;
	}

	return doom;
}


















