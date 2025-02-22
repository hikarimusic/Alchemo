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

// Update visualization style function
function updateVisualization() {
    const style = document.getElementById('styleSelect').value;
    
    // Clear and re-add the current molecule
    viewer.clear();
    
    if (currentPDBData) {
        viewer.addModel(currentPDBData, "pdb");
        // Apply selected style for protein
        viewer.setStyle({}, {[style]: {}});
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
        alert('Please select a PDB file first');
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
        
    } catch (error) {
        console.error('Error loading protein:', error);
        alert('Error loading protein file');
    }
});

// Handle drug loading
document.getElementById('loadDrugButton').addEventListener('click', async function() {
    const smileInput = document.getElementById('smileInput');
    const smileSequence = smileInput.value.trim();
    
    if (!smileSequence) {
        alert('Please enter a SMILE sequence');
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
        
    } catch (error) {
        console.error('Error loading drug:', error);
        alert('Error loading drug');
    }
});

// Handle start button
document.getElementById('startButton').addEventListener('click', async function() {
    try {
        const response = await fetch(`${API_BASE_URL}/api/get_current_state`);
        if (!response.ok) {
            throw new Error('Failed to get current state');
        }

        const data = await response.json();
        
        // Update the viewer with current state
        viewer.clear();
        if (data.pdb_data) {
            currentPDBData = data.pdb_data;
            currentDrugData = null;
        } else if (data.drug_data) {
            currentDrugData = data.drug_data;
            currentPDBData = null;
        }
        
        updateVisualization();
        
    } catch (error) {
        console.error('Error starting visualization:', error);
        alert('Error starting visualization');
    }
});

// Handle stop button
document.getElementById('stopButton').addEventListener('click', function() {
    viewer.pause();
});

// Handle style changes
document.getElementById('styleSelect').addEventListener('change', function() {
    if (currentPDBData || currentDrugData) {
        updateVisualization();
    }
});