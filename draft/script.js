// Create a 3Dmol.js viewer with a black background
var viewer = $3Dmol.createViewer("viewer", { width: "100%", height: "100%", backgroundColor: "black" });


// Function to load a local file
function loadLocalFile(file) {
    var reader = new FileReader();

    reader.onload = function(event) {
        var pdbData = event.target.result;

        viewer.clear();

        // Add the PDB data to the viewer
        viewer.addModel(pdbData, "pdb");

        // Get the selected style from the dropdown
        var styleSelect = document.getElementById("styleSelect");
        var selectedStyle = styleSelect.options[styleSelect.selectedIndex].value;

        // Apply the selected style
        if (selectedStyle === "cartoon") {
            viewer.setStyle({ cartoon: { color: 'spectrum' } });
        } else {
            // For other styles, you can define them here
            viewer.setStyle({}, { [selectedStyle]: {} });
        }

        viewer.zoomTo();
        viewer.render();
    };

    reader.readAsText(file);
}

// Event listener for the file input
document.getElementById("fileInput").addEventListener("change", function(event) {
    var fileInput = event.target;
    var selectedFile = fileInput.files[0];

    if (selectedFile) {
        loadLocalFile(selectedFile);
    }
});

// Event listener for style changes
document.getElementById("styleSelect").addEventListener("change", function() {
    var styleSelect = this;
    var selectedStyle = styleSelect.options[styleSelect.selectedIndex].value;

    // Apply the selected style
    if (selectedStyle === "cartoon") {
        viewer.setStyle({ cartoon: { color: 'spectrum' } });
    } else {
        // For other styles, you can define them here
        viewer.setStyle({}, { [selectedStyle]: {} });
    }
    viewer.render();
});