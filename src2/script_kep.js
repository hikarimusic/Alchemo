// Initialize global variables
let viewer = null;
let currentPDBData = null;
let currentDrugData = null;

// Define the API base URL
const API_BASE_URL = 'http://localhost:5000';

// Initialize the viewer when the document loads
document.addEventListener('DOMContentLoaded', function() {
    // Initialize the 3Dmol viewer
    let element = document.getElementById('viewer');
    viewer = $3Dmol.createViewer(element, {
        backgroundColor: 'black'
    });

    // Set up file input change handler
    document.getElementById('pdbFileInput').addEventListener('change', function(e) {
        const fileName = e.target.files[0] ? e.target.files[0].name : 'No file selected';
        document.getElementById('pdbFileName').value = fileName;
    });

    // Update style select options to only show supported styles
    const styleSelect = document.getElementById('styleSelect');
    styleSelect.innerHTML = `
        <option value="cartoon">Cartoon</option>
        <option value="stick">Stick</option>
        <option value="sphere">Sphere</option>
        <option value="line">Line</option>
    `;
});

// Helper function to display score
function displayScore(score) {
    const scoreDisplay = document.getElementById('scoreDisplay');
    if (score !== null) {
        scoreDisplay.textContent = `Best Score: ${score.toFixed(2)}`;
        scoreDisplay.style.display = 'block';
        
        // Auto-hide after 10 seconds
        setTimeout(() => {
            scoreDisplay.style.display = 'none';
        }, 10000);
    } else {
        scoreDisplay.style.display = 'none';
    }
}

// Helper function to show error messages
function showError(message) {
    alert(message);
    console.error(message);
}

// Update visualization style function
function updateVisualization() {
    const style = document.getElementById('styleSelect').value;
    
    // Clear and re-add the current molecule
    viewer.clear();
    
    if (currentPDBData) {
        viewer.addModel(currentPDBData, "pdb");
        
        if (style === 'cartoon') {
            // Show protein in cartoon style
            viewer.setStyle({resn: ["GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "TYR", "TRP", "SER", 
                                  "THR", "CYS", "MET", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS"]}, 
                          {cartoon: {}}); // List of standard amino acid residues
            
            // Show non-protein molecules (ligands, water, etc.) in stick style
            viewer.setStyle({hetflag: true}, {stick: {}}); // hetflag selects non-protein molecules
        } else {
            // Apply selected style to everything
            viewer.setStyle({}, {[style]: {}});
        }
    } else if (currentDrugData) {
        viewer.addModel(currentDrugData, "pdb");
        // For drugs, use stick style if cartoon is selected, otherwise use selected style
        if (style === 'cartoon') {
            viewer.setStyle({}, {stick: {}}); // Use stick style for drugs when cartoon is selected
        } else {
            viewer.setStyle({}, {[style]: {}});
        }
    }
    
    viewer.zoomTo();
    viewer.render();
}

// Handle protein file loading
document.getElementById('loadProteinButton').addEventListener('click', async function() {
    const fileInput = document.getElementById('pdbFileInput');
    const file = fileInput.files[0];
    
    if (!file) {
        showError('Please select a PDB file first');
        return;
    }

    try {
        const pdbData = await file.text();
        const response = await fetch(`${API_BASE_URL}/api/load_protein`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({ pdb_data: pdbData })
        });

        if (!response.ok) {
            throw new Error('Failed to load protein');
        }

        const result = await response.json();
        currentPDBData = result.pdb_data;
        currentDrugData = null; // Clear any existing drug data
        
        updateVisualization();
        displayScore(null); // Hide score display
        
    } catch (error) {
        showError('Error loading protein file: ' + error.message);
    }
});

// Handle drug loading
document.getElementById('loadDrugButton').addEventListener('click', async function() {
    const smileInput = document.getElementById('smileInput');
    const smileSequence = smileInput.value.trim();
    
    if (!smileSequence) {
        showError('Please enter a SMILE sequence');
        return;
    }

    try {
        const response = await fetch(`${API_BASE_URL}/api/load_drug`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({ smile_sequence: smileSequence })
        });

        if (!response.ok) {
            throw new Error('Failed to load drug');
        }

        const result = await response.json();
        currentDrugData = result.drug_data;
        currentPDBData = null; // Clear any existing protein data
        
        updateVisualization();
        displayScore(null); // Hide score display
        
    } catch (error) {
        showError('Error loading drug: ' + error.message);
    }
});

// Handle toggle button for docking
document.getElementById('toggleButton').addEventListener('click', async function() {
    const button = this;
    const currentState = button.dataset.state;

    if (currentState === 'start') {
        try {
            // First check if we have both protein and drug loaded
            const response = await fetch(`${API_BASE_URL}/api/get_current_state`);
            if (!response.ok) {
                throw new Error('Failed to get current state');
            }
            const data = await response.json();
            
            if (!data.pdb_data || !data.drug_data) {
                showError('Please load both protein and drug before starting');
                return;
            }

            // Change button state to stop
            button.textContent = 'Stop';
            button.dataset.state = 'stop';
            button.classList.add('secondary');

            // Initialize the docking process
            const initResponse = await fetch(`${API_BASE_URL}/api/initialize_docking`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                }
            });

            if (!initResponse.ok) {
                throw new Error('Failed to initialize docking');
            }

            // Start optimization loop
            let isRunning = true;
            while (isRunning && button.dataset.state === 'stop') {
                // Call optimize endpoint
                const optimizeResponse = await fetch(`${API_BASE_URL}/api/optimize_docking`, {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    body: JSON.stringify({
                        is_stopped: false
                    })
                });

                if (!optimizeResponse.ok) {
                    throw new Error('Failed to optimize conformations');
                }

                const result = await optimizeResponse.json();
                
                // Update visualization
                viewer.clear();
                
                // Add protein structure
                viewer.addModel(result.protein_pdb, "pdb");
                
                // Get current style
                const style = document.getElementById('styleSelect').value;
                
                if (style === 'cartoon') {
                    // Show protein in cartoon style
                    viewer.setStyle({resn: ["GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "TYR", "TRP", "SER",
                                         "THR", "CYS", "MET", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS"]},
                                {cartoon: {}});
                    
                    // Any non-protein parts in the protein model should be in stick style
                    viewer.setStyle({hetflag: true}, {stick: {}});
                } else {
                    // Apply selected style to all atoms
                    viewer.setStyle({}, {[style]: {}});
                }
                
                // Add each drug conformer as a separate model with different colors
                result.conformer_pdbs.forEach((conformerPdb, index) => {
                    const model = viewer.addModel(conformerPdb, "pdb");
                    const t = index / Math.max(1, result.conformer_pdbs.length - 1); // normalized position in spectrum
                    const r = Math.floor(255 * (1 - t)); // red component decreases
                    const b = Math.floor(255 * t);       // blue component increases
                    const color = `0x${r.toString(16).padStart(2,'0')}00${b.toString(16).padStart(2,'0')}`; // format as hex
                    
                    if (style === 'cartoon') {
                        viewer.setStyle({model: model}, {stick: {color: color}});
                    } else {
                        viewer.setStyle({model: model}, {[style]: {color: color}});
                    }
                });
                
                viewer.zoomTo();
                viewer.render();

                // During optimization, don't display score
                displayScore(null);

                // Check if optimization is complete
                if (result.is_complete) {
                    isRunning = false;
                }

                // Add a small delay to prevent overwhelming the server
                await new Promise(resolve => setTimeout(resolve, 100));
            }

            // If we broke out of the loop due to stop button, call optimize once more with is_stopped=true
            if (button.dataset.state === 'start') {
                // User pressed stop - retrieve only the best conformer
                const finalResponse = await fetch(`${API_BASE_URL}/api/optimize_docking`, {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    body: JSON.stringify({
                        is_stopped: true
                    })
                });

                if (finalResponse.ok) {
                    const finalResult = await finalResponse.json();
                    
                    // Update visualization with only the best conformer
                    viewer.clear();
                    
                    // Add protein structure
                    viewer.addModel(finalResult.protein_pdb, "pdb");
                    
                    // Get current style
                    const style = document.getElementById('styleSelect').value;
                    
                    if (style === 'cartoon') {
                        // Show protein in cartoon style
                        viewer.setStyle({resn: ["GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "TYR", "TRP", "SER",
                                             "THR", "CYS", "MET", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS"]},
                                    {cartoon: {}});
                        
                        // Any non-protein parts in the protein model should be in stick style
                        viewer.setStyle({hetflag: true}, {stick: {}});
                    } else {
                        // Apply selected style to all atoms
                        viewer.setStyle({}, {[style]: {}});
                    }
                    
                    // Add only the best conformer in a distinct color
                    if (finalResult.conformer_pdbs.length > 0) {
                        const model = viewer.addModel(finalResult.conformer_pdbs[0], "pdb");
                        
                        // Use bright green for the best conformer
                        const color = "0x00FF00"; // Bright green
                        
                        if (style === 'cartoon') {
                            viewer.setStyle({model: model}, {stick: {color: color}});
                        } else {
                            viewer.setStyle({model: model}, {[style]: {color: color}});
                        }
                        
                        // Display the score
                        if (finalResult.scores && finalResult.scores.length > 0) {
                            displayScore(finalResult.scores[0]);
                        }
                    }
                    
                    viewer.zoomTo();
                    viewer.render();
                }
            }

            // Reset button state
            button.textContent = 'Start';
            button.dataset.state = 'start';
            button.classList.remove('secondary');

        } catch (error) {
            showError('Error during docking: ' + error.message);
            
            // Reset button state on error
            button.textContent = 'Start';
            button.dataset.state = 'start';
            button.classList.remove('secondary');
        }
    } else {
        // User clicked stop - button state will be checked in while loop
        button.textContent = 'Start';
        button.dataset.state = 'start';
        button.classList.remove('secondary');
    }
});

// Handle style changes
document.getElementById('styleSelect').addEventListener('change', function() {
    if (currentPDBData || currentDrugData) {
        updateVisualization();
    }
});