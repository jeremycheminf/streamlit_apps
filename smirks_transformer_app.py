#!/usr/bin/env python3
"""
Streamlit app for SMIRKS transformation of SMILES molecules.

This app allows users to input a SMILES string and a SMIRKS transformation,
then displays the original and transformed molecules as images.

Requirements:
- RDKit (available in cheminf_utils environment)
- Streamlit
"""

import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdChemReactions
import base64
from io import BytesIO
import re
import tempfile
import os

# Set page configuration
st.set_page_config(
    page_title="SMIRKS Transformer",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS for better styling
st.markdown(
    """
<style>
    .main-header {
        font-size: 3rem;
        font-weight: bold;
        text-align: center;
        margin-bottom: 2rem;
        color: #1f77b4;
    }
    .sub-header {
        font-size: 1.5rem;
        font-weight: bold;
        margin-top: 2rem;
        margin-bottom: 1rem;
        color: #333;
    }
    .info-box {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin-bottom: 2rem;
    }
    .error-box {
        background-color: #ffebee;
        border: 1px solid #f44336;
        color: #d32f2f;
        padding: 1rem;
        border-radius: 0.5rem;
        margin-bottom: 1rem;
    }
    .success-box {
        background-color: #e8f5e9;
        border: 1px solid #4caf50;
        color: #2e7d32;
        padding: 1rem;
        border-radius: 0.5rem;
        margin-bottom: 1rem;
    }
</style>
""",
    unsafe_allow_html=True,
)


def smiles_to_mol(smiles):
    """Convert SMILES string to RDKit molecule."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, "Invalid SMILES string"
        return mol, None
    except Exception as e:
        return None, f"Error parsing SMILES: {str(e)}"


def smirks_to_rxn(smirks):
    """Convert SMIRKS string to RDKit reaction."""
    try:
        # RDKit expects SMARTS reaction format: reactants > agents > products
        # Convert common SMIRKS format (with >>) to SMARTS format (with >)
        if ">>" in smirks and ">" not in smirks.replace(">>", ""):
            # Convert >> format to > format by adding empty agent
            smirks = smirks.replace(">>", "> >")

        rxn = rdChemReactions.ReactionFromSmarts(smirks)
        if rxn is None:
            return None, "Invalid reaction SMARTS string"
        return rxn, None
    except Exception as e:
        return None, f"Error parsing reaction: {str(e)}"


def apply_smirks_to_molecule(mol, rxn):
    """Apply SMIRKS transformation to a molecule."""
    try:
        products = rxn.RunReactants([mol])
        if not products:
            return None, "No products generated - SMIRKS may not match the molecule"

        # Return the first product
        product_mol = products[0][0]
        if product_mol is None:
            return None, "Failed to generate product molecule"

        return product_mol, None
    except Exception as e:
        return None, f"Error applying SMIRKS: {str(e)}"


def mol_to_image_base64(mol, size=(400, 400)):
    """Convert RDKit molecule to base64 encoded PNG image."""
    try:
        if mol is None:
            return None

        # Generate 2D coordinates if not present
        if not mol.GetNumConformers():
            AllChem.Compute2DCoords(mol)

        # Create drawing
        drawer = rdMolDraw2D.MolDraw2DCairo(*size)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        # Convert to base64
        png_data = drawer.GetDrawingText()
        b64_string = base64.b64encode(png_data).decode("utf-8")

        return b64_string
    except Exception as e:
        st.error(f"Error generating image: {str(e)}")
        return None


def validate_smirks_format(smirks):
    """Basic validation of reaction format."""
    # Check for required arrow (either >> or >)
    if ">>" not in smirks and ">" not in smirks:
        return (
            False,
            "Reaction must contain '>>' or '>' to separate reactants from products",
        )

    # Check for balanced parentheses (basic check)
    if smirks.count("(") != smirks.count(")"):
        return False, "Unbalanced parentheses in reaction"

    return True, ""


def rxn_file_to_smirks(rxn_file_content):
    """Convert RXN file content to SMIRKS string."""
    try:
        # Write the RXN content to a temporary file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".rxn", delete=False
        ) as temp_file:
            temp_file.write(rxn_file_content)
            temp_file_path = temp_file.name

        try:
            # Read the reaction from the RXN file
            rxn = AllChem.ReactionFromRxnFile(temp_file_path, sanitize=False)
            if rxn is None:
                return None, "Failed to parse RXN file"

            # Convert to SMARTS
            smirks = AllChem.ReactionToSmarts(rxn)
            return smirks, None
        finally:
            # Clean up the temporary file
            os.unlink(temp_file_path)

    except Exception as e:
        return None, f"Error processing RXN file: {str(e)}"


def process_rxn_file(uploaded_file):
    """Process uploaded RXN file and return SMIRKS string."""
    try:
        # Read the file content
        content = uploaded_file.getvalue().decode("utf-8")

        # Convert to SMIRKS
        smirks, error = rxn_file_to_smirks(content)
        if error:
            return None, error

        return smirks, None
    except Exception as e:
        return None, f"Error reading file: {str(e)}"


def main():
    # Header
    st.markdown(
        '<div class="main-header">🧪 SMIRKS Transformer</div>', unsafe_allow_html=True
    )

    # Sidebar with instructions
    with st.sidebar:
        st.header("Instructions")
        st.markdown(
            """
        1. **Enter a SMILES string** - The molecule to transform
        2. **Enter a reaction SMARTS string** - The transformation rule
        3. **OR upload an RXN file** - Automatically convert to SMIRKS
        4. **Click "Transform"** - See the result
        
        **Examples:**
        - SMILES: `CC(=O)OC1=CC=CC=C1C(=O)O` (Aspirin)
        - Reaction: `[C:1][O:2]>>[C:1][N:2]` (Replace O with N)
        - Reaction: `[c:1][O:2]>>[c:1][Cl:2]` (Replace O with Cl on aromatic ring)
        """
        )

        st.markdown("---")
        st.header("RXN File Support")
        st.markdown(
            """
        Upload an RXN file to automatically extract the reaction and populate the SMIRKS field.
        
        **Supported formats:**
        - RXN files (Reaction files)
        - Automatic conversion to SMIRKS format
        """
        )

        st.markdown("---")
        st.header("About")
        st.markdown(
            """
        This app uses RDKit to:
        - Parse SMILES and reaction SMARTS strings
        - Apply chemical transformations
        - Generate 2D molecular depictions
        """
        )

    # Main content
    col1, col2 = st.columns([1, 1])

    with col1:
        st.markdown(
            '<div class="sub-header">Input Molecule</div>', unsafe_allow_html=True
        )

        # SMILES input
        smiles_input = st.text_area(
            "SMILES String:",
            value="CC(=O)OC1=CC=CC=C1C(=O)O",
            height=100,
            help="Enter a valid SMILES string for the molecule you want to transform",
        )

        st.markdown(
            '<div class="sub-header">Reaction Input</div>', unsafe_allow_html=True
        )

        # RXN file upload
        rxn_file = st.file_uploader(
            "Upload RXN File (optional)",
            type=["rxn"],
            help="Upload an RXN file to automatically extract the reaction and populate the SMIRKS field",
        )

        # SMIRKS input
        smirks_input = st.text_area(
            "SMIRKS Transformation:",
            value="[C:1][O:2]>>[C:1][N:2]",
            height=100,
            help="Enter a valid SMIRKS string for the transformation",
        )

        # Process RXN file if uploaded
        if rxn_file is not None:
            with st.spinner("Processing RXN file..."):
                smirks_from_file, error = process_rxn_file(rxn_file)
                if error:
                    st.error(f"Error processing RXN file: {error}")
                else:
                    # Update the SMIRKS input with the extracted SMIRKS
                    st.session_state.smirks_input = smirks_from_file
                    st.success(f"Successfully extracted SMIRKS from RXN file")
                    st.code(smirks_from_file, language="text")
                    # Auto-populate the SMIRKS field
                    smirks_input = smirks_from_file

        # Transform button
        transform_button = st.button("🧪 Transform Molecule", type="primary")

    with col2:
        st.markdown('<div class="sub-header">Results</div>', unsafe_allow_html=True)

        if transform_button:
            # Validate inputs
            if not smiles_input.strip():
                st.markdown(
                    '<div class="error-box">❌ Please enter a SMILES string</div>',
                    unsafe_allow_html=True,
                )
                return

            if not smirks_input.strip():
                st.markdown(
                    '<div class="error-box">❌ Please enter a SMIRKS string</div>',
                    unsafe_allow_html=True,
                )
                return

            # Validate SMIRKS format
            is_valid, error_msg = validate_smirks_format(smirks_input)
            if not is_valid:
                st.markdown(
                    f'<div class="error-box">❌ {error_msg}</div>',
                    unsafe_allow_html=True,
                )
                return

            # Process SMILES
            with st.spinner("Processing SMILES..."):
                mol, error = smiles_to_mol(smiles_input.strip())
                if error:
                    st.markdown(
                        f'<div class="error-box">❌ SMILES Error: {error}</div>',
                        unsafe_allow_html=True,
                    )
                    return

            # Process SMIRKS
            with st.spinner("Processing SMIRKS..."):
                rxn, error = smirks_to_rxn(smirks_input.strip())
                if error:
                    st.markdown(
                        f'<div class="error-box">❌ SMIRKS Error: {error}</div>',
                        unsafe_allow_html=True,
                    )
                    return

            # Apply transformation
            with st.spinner("Applying transformation..."):
                product_mol, error = apply_smirks_to_molecule(mol, rxn)
                if error:
                    st.markdown(
                        f'<div class="error-box">❌ Transformation Error: {error}</div>',
                        unsafe_allow_html=True,
                    )
                    return

            # Generate images
            with st.spinner("Generating molecular images..."):
                original_img = mol_to_image_base64(mol)
                product_img = mol_to_image_base64(product_mol)

            if original_img and product_img:
                # Display results
                st.markdown(
                    '<div class="success-box">✅ Transformation successful!</div>',
                    unsafe_allow_html=True,
                )

                # Show SMILES strings
                original_smiles = Chem.MolToSmiles(mol, canonical=True)
                product_smiles = Chem.MolToSmiles(product_mol, canonical=True)

                st.markdown("### Original Molecule")
                st.code(original_smiles, language="text")

                st.markdown("### Transformed Molecule")
                st.code(product_smiles, language="text")

                # Display images side by side
                img_col1, img_col2 = st.columns([1, 1])

                with img_col1:
                    st.markdown("#### Before")
                    st.image(
                        f"data:image/png;base64,{original_img}", use_column_width=True
                    )

                with img_col2:
                    st.markdown("#### After")
                    st.image(
                        f"data:image/png;base64,{product_img}", use_column_width=True
                    )

                # Download options
                st.markdown("### Download Results")
                col_download1, col_download2 = st.columns([1, 1])

                with col_download1:
                    st.download_button(
                        label="Download Original SMILES",
                        data=original_smiles,
                        file_name="original_molecule.smiles",
                        mime="text/plain",
                    )

                with col_download2:
                    st.download_button(
                        label="Download Transformed SMILES",
                        data=product_smiles,
                        file_name="transformed_molecule.smiles",
                        mime="text/plain",
                    )
            else:
                st.markdown(
                    '<div class="error-box">❌ Failed to generate molecular images</div>',
                    unsafe_allow_html=True,
                )
        else:
            # Default display when no transformation has been performed
            st.markdown(
                '<div class="info-box">💡 Enter SMILES and SMIRKS strings above, then click "Transform" to see the results.</div>',
                unsafe_allow_html=True,
            )

            # Show example transformation
            if st.button("Show Example"):
                # Use the default values to show an example
                mol, _ = smiles_to_mol("CC(=O)OC1=CC=CC=C1C(=O)O")
                rxn, _ = smirks_to_rxn("[C:1][O:2]>>[C:1][N:2]")
                product_mol, _ = apply_smirks_to_molecule(mol, rxn)

                original_img = mol_to_image_base64(mol)
                product_img = mol_to_image_base64(product_mol)

                if original_img and product_img:
                    st.markdown("### Example Transformation")

                    img_col1, img_col2 = st.columns([1, 1])

                    with img_col1:
                        st.markdown("#### Before")
                        st.image(
                            f"data:image/png;base64,{original_img}",
                            use_column_width=True,
                        )

                    with img_col2:
                        st.markdown("#### After")
                        st.image(
                            f"data:image/png;base64,{product_img}",
                            use_column_width=True,
                        )


if __name__ == "__main__":
    main()
