// Create a 3Dmol.js viewer with a black background
var viewer = $3Dmol.createViewer("viewer", { width: "100%", height: "100%", backgroundColor: "black" });

// Function to load a local file
function loadLocalFile(file) {
    var reader = new FileReader();

    reader.onload = function(event) {
        var pdbData = event.target.result;
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
    $3Dmol.download("pdb:" + pdbId, viewer, {}, function(m) {
        if (m) {
            applySelectedStyle('protein');
            viewer.zoomTo();
            viewer.render();
        } else {
            alert("Error loading PDB ID.");
        }
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
// function applySelectedStyle() {
//     var styleSelect = document.getElementById("styleSelect");
//     var selectedStyle = styleSelect.options[styleSelect.selectedIndex].value;
//     var styleOptions = { stick: {}, sphere: { scale: 0.3 }, cartoon: { color: 'spectrum' } };

//     viewer.setStyle({}, styleOptions[selectedStyle] || {});
// }

function applySelectedStyle(moleculeType) {
    var styleSelect = document.getElementById("styleSelect");
    var selectedStyle = styleSelect.options[styleSelect.selectedIndex].value;

    // Apply the selected style
    if (moleculeType === "protein" && selectedStyle === "cartoon") {
        viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
    } else {
        // For other styles, you can define them here
        viewer.setStyle({}, { [selectedStyle]: {} });
    }
}

// Event listeners
document.getElementById("loadFileButton").addEventListener("click", function() {
    var fileInput = document.getElementById("fileInput");
    var file = fileInput.files[0];
    if (file) {
        loadLocalFile(file);
    }
});

// Event listener for the file input
// document.getElementById("fileInput").addEventListener("change", function(event) {
//     var fileInput = event.target;
//     var selectedFile = fileInput.files[0];

//     if (selectedFile) {
//         loadLocalFile(selectedFile);
//     }
// });

// Event listener for style changes

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

document.getElementById("clearViewButton").addEventListener("click", function() {
    viewer.clear();
    viewer.render();
});
