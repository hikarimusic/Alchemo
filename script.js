// Create a 3Dmol.js viewer with a black background
var viewer = $3Dmol.createViewer("viewer", { width: "100%", height: "100%", backgroundColor: "black" });

// Save th pdbData for simulation
var pdbData = ""

// Function to load a local file
function loadLocalFile(file) {
    var reader = new FileReader();

    reader.onload = function(event) {
        pdbData = event.target.result;
        viewer.clear();
        viewer.addModel(pdbData, "pdb");
        applySelectedStyle('protein');
        viewer.zoomTo();
        viewer.render();
    };

    reader.readAsText(file);
}

// Function to load a structure from PDB ID
function loadFromPdbId(pdbId) {
    viewer.clear();

    var pdbUrl = `https://files.rcsb.org/download/${pdbId}.pdb`;
    fetch(pdbUrl)
        .then(response => {
            if (!response.ok) {
                throw new Error('Network response was not ok');
            }
            return response.text();
        })
        .then(data => {
            pdbData = data;
            viewer.addModel(pdbData, "pdb");
            applySelectedStyle('protein');
            viewer.zoomTo();
            viewer.render();
        })
        .catch(error => {
            console.error('There has been a problem with your fetch operation:', error);
            alert("Error loading PDB ID.");
        });
}

// Function to load a structure from PubChem ID
function loadFromPubchemId(pubchemId) {
    viewer.clear();
    $3Dmol.download("cid:" + pubchemId, viewer, {}, function(m) {
        if (m) {
            applySelectedStyle('compound');
            viewer.zoomTo();
            viewer.render();
        } else {
            alert("Error loading PubChem ID.");
        }
    });
}

// Function to apply the selected style
function applySelectedStyle(moleculeType) {
    var styleSelect = document.getElementById("styleSelect");
    var selectedStyle = styleSelect.options[styleSelect.selectedIndex].value;

    if (moleculeType === "protein" && selectedStyle === "cartoon") {
        viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
    } else if (moleculeType === "compound" && selectedStyle === "cartoon") {
        viewer.setStyle({}, { stick: {} });
    } else {
        viewer.setStyle({}, { [selectedStyle]: {} });
    }
}


// Function to animate and stop the protein
var proteinFrames = [];
var currentFrameIndex = 0;
var continueAnimation = true;

function animateProtein() {
    fetch('http://localhost:5000/animateProtein', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({ pdbData: pdbData })
    })
    .then(response => response.json())
    .then(data => {
        proteinFrames = data.frames;
        animateProteinFrames();
    })
    .catch(error => console.error('Error sending protein data to backend:', error));
    continueAnimation = true;
}

function animateProteinFrames() {
    if (continueAnimation && currentFrameIndex < proteinFrames.length) {
        viewer.clear();
        viewer.addModel(proteinFrames[currentFrameIndex], "pdb");
        applySelectedStyle('protein');
        viewer.render();

        currentFrameIndex++;
        setTimeout(animateProteinFrames, 100);
    } else {
        currentFrameIndex = 0;
        animateProteinFrames();
    }
}

function stopAnimation() {
    continueAnimation = false;
}

// Event listeners
document.getElementById("loadFileButton").addEventListener("click", function() {
    var fileInput = document.getElementById("fileInput");
    var file = fileInput.files[0];
    if (file) {
        loadLocalFile(file);
    }
});

document.getElementById("loadPdbIdButton").addEventListener("click", function() {
    var pdbId = document.getElementById("pdbIdInput").value;
    loadFromPdbId(pdbId);
});

document.getElementById("loadPubchemIdButton").addEventListener("click", function() {
    var pubchemId = document.getElementById("pubchemIdInput").value;
    loadFromPubchemId(pubchemId);
});

document.getElementById("styleSelect").addEventListener("change", function() {
    var styleSelect = this;
    var selectedStyle = styleSelect.options[styleSelect.selectedIndex].value;

    applySelectedStyle('protein');
    viewer.render();
});

document.getElementById("animateProteinButton").addEventListener("click", function() {
    animateProtein();
});

document.getElementById("stopAnimationButton").addEventListener("click", function() {
    stopAnimation();
});
