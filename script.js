// 3Dmol Viewer 
var viewer = $3Dmol.createViewer("viewer", { width: "100%", height: "100%", backgroundColor: "black" });

// Load Protein
function loadProtein() {
    var proteinSequence = document.getElementById("sequenceInput").value;
    fetch('http://localhost:5000/loadProtein', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({ sequence: proteinSequence })
    })
    .then(response => response.json())
    .then(() => showProtein())
    .catch(error => console.error('Error:', error));
}

// Load PdbId
function loadPdbId() {
    var pdbId = document.getElementById("idInput").value;
    var pdbUrl = `https://files.rcsb.org/download/${pdbId}.pdb`;
    fetch(pdbUrl)
    .then(response => {
        if (!response.ok) {
            throw new Error('Network response was not ok ' + response.statusText);
        }
        return response.text();
    })
    .then(pdbData => {
        fetch('http://localhost:5000/loadPdbId', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ pdbData: pdbData })
        })
        .then(response => response.json())
        .then(() => showProtein())
        .catch(error => console.error('Error:', error));
    })
    .catch(error => console.error('Error fetching PDB file:', error));
}

// Apply Style
function applyStyle(moleculeType) {
    var styleSelect = document.getElementById("styleSelect");
    var selectedStyle = styleSelect.options[styleSelect.selectedIndex].value;

    if (moleculeType === "protein" && selectedStyle === "cartoon") {
        viewer.setStyle({ chain: 'A' }, { cartoon: { color: 'spectrum' } });
        viewer.setStyle({ chain: 'B' }, { sphere: { radius: 0.3, color: 'aqua' } });
    } else if (moleculeType === "compound" && selectedStyle === "cartoon") {
        viewer.setStyle({ chain: 'A' }, { stick: {} });
        viewer.setStyle({ chain: 'B' }, { sphere: { radius: 0.3, color: 'aqua' } });
    } else {
        viewer.setStyle({ chain: 'A' }, { [selectedStyle]: {} });
        viewer.setStyle({ chain: 'B' }, { sphere: { radius: 0.3, color: 'aqua' } });
    }
}

// Show Protein
function showProtein() {
    fetch('http://localhost:5000/showProtein')
    .then(response => response.json())
    .then(data => {
        viewer.clear();
        viewer.addModel(data.pdbData, "pdb");
        applyStyle('protein');
        viewer.zoomTo();
        viewer.render();
    })
    .catch(error => console.error('Error fetching protein structure:', error));
}

// Simulate Portein
function simulateProtein() {
    fetch('http://localhost:5000/simulateProtein')
    .then(response => response.json())
    .then(data => {
        viewer.clear();
        viewer.addModel(data.pdbData, "pdb");
        applyStyle('protein');
        // viewer.zoomTo();
        viewer.render();
        if (isSimulationRunning) {
            setTimeout(simulateProtein, 1);
        }
    })
    .catch(error => {
        console.error('Error simulating protein structure:', error);
        isSimulationRunning = false;
    });
}

// Event Listeners
document.getElementById("loadProteinButton").addEventListener("click", loadProtein);

document.getElementById("loadPdbIdButton").addEventListener("click", loadPdbId)

document.getElementById("styleSelect").addEventListener("change", function() {
    var styleSelect = this;
    var selectedStyle = styleSelect.options[styleSelect.selectedIndex].value;
    applyStyle('protein');
    viewer.render();
});

document.getElementById("simulateButton").addEventListener("click", function() {
    isSimulationRunning = true;
    showProtein()
    simulateProtein();
});

document.getElementById("stopButton").addEventListener("click", function() {
    isSimulationRunning = false;
});



























// function simulateProtein() {
//     fetch('http://localhost:5000/simulateProtein')
//     .then(response => response.json())
//     .then(data => {
//         // pdbData = data.pdbData;
//         viewer.clear();
//         viewer.addModel(data.pdbData, "pdb");
//         applyStyle('protein');
//         viewer.zoomTo();
//         viewer.render();
//     })
//     .catch(error => console.error('Error fetching protein structure:', error));
// }


// function animateProteinFrames() {
//     if (continueAnimation && currentFrameIndex < proteinFrames.length) {
//         viewer.clear();
//         viewer.addModel(proteinFrames[currentFrameIndex], "pdb");
//         applyStyle('protein');
//         viewer.render();

//         currentFrameIndex++;
//         setTimeout(animateProteinFrames, 100);
//     } else {
//         currentFrameIndex = 0;
//         animateProteinFrames();
//     }
// }

// function stopProtein() {
//     continueAnimation = false;
// }




// var proteinFrames = [];
// var currentFrameIndex = 0;
// var continueAnimation = true;

// function simulateProtein() {
//     fetch('http://localhost:5000/simulateProtein', {
//         method: 'POST',
//         headers: {
//             'Content-Type': 'application/json'
//         },
//         body: JSON.stringify({ pdbData: pdbData })
//     })
//     .then(response => response.json())
//     .then(data => {
//         proteinFrames = data.frames;
//         animateProteinFrames();
//     })
//     .catch(error => console.error('Error sending protein data to backend:', error));
//     continueAnimation = true;
// }

// function animateProteinFrames() {
//     if (continueAnimation && currentFrameIndex < proteinFrames.length) {
//         viewer.clear();
//         viewer.addModel(proteinFrames[currentFrameIndex], "pdb");
//         applyStyle('protein');
//         viewer.render();

//         currentFrameIndex++;
//         setTimeout(animateProteinFrames, 100);
//     } else {
//         currentFrameIndex = 0;
//         animateProteinFrames();
//     }
// }

// function stopProtein() {
//     continueAnimation = false;
// }


// Function to load protein sequence
// function loadProteinSequence() {
//     // var proteinSequence = document.getElementById("proteinInput").value;
//     var proteinSequence = document.getElementById("sequenceInput").value;
//     fetch('http://localhost:5000/generateProteinStructure', {
//         method: 'POST',
//         headers: {
//             'Content-Type': 'application/json'
//         },
//         body: JSON.stringify({ sequence: proteinSequence })
//     })
//     .then(response => response.json())
//     .then(data => {
//         pdbData = data.pdbData;
//         viewer.clear();
//         viewer.addModel(pdbData, "pdb");
//         applySelectedStyle('protein');
//         viewer.zoomTo();
//         viewer.render();
//     })
//     .catch(error => console.error('Error:', error));
// }

// Function to load molecule sequence
// function loadMoleculeSequence() {
//     // var moleculeSequence = document.getElementById("moleculeInput").value;
//     var moleculeSequence = document.getElementById("sequenceInput").value;
//     fetch('http://localhost:5000/generateMoleculeStructure', {
//         method: 'POST',
//         headers: {
//             'Content-Type': 'application/json'
//         },
//         body: JSON.stringify({ sequence: moleculeSequence })
//     })
//     .then(response => response.json())
//     .then(data => {
//         pdbData = data.pdbData;
//         viewer.clear();
//         viewer.addModel(pdbData, "pdb");
//         applyStyle('compound');
//         viewer.zoomTo();
//         viewer.render();
//     })
//     .catch(error => console.error('Error:', error));
// }

// Function to load a structure from PDB ID
// function loadFromPdbId(pdbId) {
//     viewer.clear();

//     var pdbUrl = `https://files.rcsb.org/download/${pdbId}.pdb`;
//     fetch(pdbUrl)
//         .then(response => {
//             if (!response.ok) {
//                 throw new Error('Network response was not ok');
//             }
//             return response.text();
//         })
//         .then(data => {
//             pdbData = data;
//             viewer.addModel(pdbData, "pdb");
//             applyStyle('protein');
//             viewer.zoomTo();
//             viewer.render();
//         })
//         .catch(error => {
//             console.error('There has been a problem with your fetch operation:', error);
//             alert("Error loading PDB ID.");
//         });
// }

// Function to load a structure from PubChem ID
// function loadFromPubchemId(pubchemId) {
//     viewer.clear();
//     $3Dmol.download("cid:" + pubchemId, viewer, {}, function(m) {
//         if (m) {
//             applyStyle('compound');
//             viewer.zoomTo();
//             viewer.render();
//         } else {
//             alert("Error loading PubChem ID.");
//         }
//     });
// }

// Function to load a local file
// function loadLocalFile(file) {
//     var reader = new FileReader();

//     reader.onload = function(event) {
//         pdbData = event.target.result;
//         viewer.clear();
//         viewer.addModel(pdbData, "pdb");
//         applyStyle('protein');
//         viewer.zoomTo();
//         viewer.render();
//     };

//     reader.readAsText(file);
// }

// Event listeners

// document.getElementById("loadProteinButton").addEventListener("click", loadProteinSequence);

// document.getElementById("loadMoleculeButton").addEventListener("click", loadMoleculeSequence);

// document.getElementById("loadPdbIdButton").addEventListener("click", function() {
//     // var pdbId = document.getElementById("pdbIdInput").value;
//     var pdbId = document.getElementById("idInput").value;
//     loadFromPdbId(pdbId);
// });

// document.getElementById("loadPubchemIdButton").addEventListener("click", function() {
//     // var pubchemId = document.getElementById("pubchemIdInput").value;
//     var pubchemId = document.getElementById("idInput").value;
//     loadFromPubchemId(pubchemId);
// });

// document.getElementById("loadFileButton").addEventListener("click", function() {
//     var fileInput = document.getElementById("fileInput");
//     var file = fileInput.files[0];
//     if (file) {
//         loadLocalFile(file);
//     }
// });

// document.getElementById("styleSelect").addEventListener("change", function() {
//     var styleSelect = this;
//     var selectedStyle = styleSelect.options[styleSelect.selectedIndex].value;

//     applySelectedStyle('protein');
//     viewer.render();
// });

// document.getElementById("animateProteinButton").addEventListener("click", function() {
//     animateProtein();
// });

// document.getElementById("stopAnimationButton").addEventListener("click", function() {
//     stopAnimation();
// });
