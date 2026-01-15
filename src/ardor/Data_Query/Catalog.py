# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 2026

@author: natha

Catalog class for easy access to catalog files and their parameters.
"""

import pandas as pd
import os
import glob


class Catalog:
    """
    A class to manage and access catalog files in the catalog directory.
    
    This class provides convenient methods to load, list, and query catalog files,
    making it easy to retrieve parameters from available catalogs.
    
    Attributes
    ----------
    catalog_dir : str
        Path to the catalog directory.
    catalogs : dict
        Dictionary mapping catalog names to their loaded DataFrames.
    """
    
    def __init__(self, catalog_dir="/ugrad/whitsett.n/ardor/catalogs"):
        """
        Initialize the Catalog class.
        
        Parameters
        ----------
        catalog_dir : str, optional
            Path to the catalog directory. Default is "/ugrad/whitsett.n/ardor/catalogs".
        """
        self.catalog_dir = catalog_dir
        self.catalogs = {}
        
        # Check if directory exists
        if not os.path.exists(self.catalog_dir):
            raise ValueError(f"Catalog directory does not exist: {self.catalog_dir}")
        
        print(f"✅ Catalog initialized. Directory: {self.catalog_dir}")
    
    def list_available_catalogs(self):
        """
        List all available catalog files in the catalog directory.
        
        Returns
        -------
        list
            List of catalog filenames (including extensions).
        """
        catalog_files = glob.glob(os.path.join(self.catalog_dir, "*.csv"))
        catalog_names = [os.path.basename(f) for f in catalog_files]
        
        print(f"📋 Available catalogs: {len(catalog_names)}")
        for name in sorted(catalog_names):
            print(f"   - {name}")
        
        return sorted(catalog_names)
    
    def load_catalog(self, catalog_name):
        """
        Load a catalog file into memory.
        
        Parameters
        ----------
        catalog_name : str
            Name of the catalog file (e.g., 'planet_params.csv').
        
        Returns
        -------
        pd.DataFrame
            The loaded catalog as a DataFrame.
        """
        catalog_path = os.path.join(self.catalog_dir, catalog_name)
        
        if not os.path.exists(catalog_path):
            raise FileNotFoundError(f"Catalog file not found: {catalog_path}")
        
        try:
            df = pd.read_csv(catalog_path)
            self.catalogs[catalog_name] = df
            print(f"✅ Loaded catalog '{catalog_name}' with {len(df)} rows")
            return df
        except Exception as e:
            print(f"❌ Failed to load catalog '{catalog_name}'. Error: {e}")
            return None
    
    def get_catalog(self, catalog_name):
        """
        Get a catalog, loading it if necessary.
        
        Parameters
        ----------
        catalog_name : str
            Name of the catalog file.
        
        Returns
        -------
        pd.DataFrame
            The catalog as a DataFrame.
        """
        if catalog_name not in self.catalogs:
            return self.load_catalog(catalog_name)
        return self.catalogs[catalog_name]
    
    def query_catalog(self, catalog_name, **filters):
        """
        Query a catalog with specified filter conditions.
        
        Parameters
        ----------
        catalog_name : str
            Name of the catalog file.
        **filters : keyword arguments
            Filtering conditions where key is column name and value is the filter value.
            Example: hostname='AU Mic', planet_letter='b'
        
        Returns
        -------
        pd.DataFrame
            Filtered catalog rows matching all conditions.
        """
        df = self.get_catalog(catalog_name)
        
        if df is None:
            return None
        
        result = df.copy()
        
        for column, value in filters.items():
            if column not in result.columns:
                print(f"⚠️  Column '{column}' not found in catalog")
                continue
            result = result[result[column] == value]
        
        print(f"✅ Query returned {len(result)} matching rows")
        return result
    
    def get_host_data(self, catalog_name, hostname):
        """
        Get all data for a specific host from a catalog.
        
        Parameters
        ----------
        catalog_name : str
            Name of the catalog file.
        hostname : str
            Name of the host to retrieve.
        
        Returns
        -------
        pd.DataFrame
            All rows matching the hostname.
        """
        return self.query_catalog(catalog_name, hostname=hostname)
    
    def get_planet_data(self, catalog_name, hostname, planet_letter):
        """
        Get data for a specific planet around a host.
        
        Parameters
        ----------
        catalog_name : str
            Name of the catalog file.
        hostname : str
            Name of the host star.
        planet_letter : str
            Letter designation of the planet (e.g., 'b', 'c').
        
        Returns
        -------
        pd.DataFrame
            Row(s) matching the hostname and planet letter.
        """
        return self.query_catalog(catalog_name, hostname=hostname, planet_letter=planet_letter)
    
    def get_column_values(self, catalog_name, column_name, **filters):
        """
        Get specific column values from a catalog with optional filters.
        
        Parameters
        ----------
        catalog_name : str
            Name of the catalog file.
        column_name : str
            Name of the column to retrieve values from.
        **filters : keyword arguments
            Optional filtering conditions.
        
        Returns
        -------
        list
            List of values from the specified column.
        """
        if filters:
            result = self.query_catalog(catalog_name, **filters)
        else:
            result = self.get_catalog(catalog_name)
        
        if result is None or len(result) == 0:
            return []
        
        if column_name not in result.columns:
            print(f"⚠️  Column '{column_name}' not found in catalog")
            return []
        
        return result[column_name].tolist()
    
    def get_unique_hosts(self, catalog_name):
        """
        Get list of unique hosts in a catalog.
        
        Parameters
        ----------
        catalog_name : str
            Name of the catalog file.
        
        Returns
        -------
        list
            List of unique hostnames.
        """
        df = self.get_catalog(catalog_name)
        
        if df is None or 'hostname' not in df.columns:
            print("⚠️  'hostname' column not found in catalog")
            return []
        
        unique_hosts = df['hostname'].unique().tolist()
        print(f"✅ Found {len(unique_hosts)} unique hosts: {unique_hosts}")
        return unique_hosts
    
    def get_catalog_info(self, catalog_name):
        """
        Get information about a catalog (column names, size, etc.).
        
        Parameters
        ----------
        catalog_name : str
            Name of the catalog file.
        
        Returns
        -------
        dict
            Dictionary with catalog information.
        """
        df = self.get_catalog(catalog_name)
        
        if df is None:
            return None
        
        info = {
            'filename': catalog_name,
            'filepath': os.path.join(self.catalog_dir, catalog_name),
            'num_rows': len(df),
            'num_columns': len(df.columns),
            'columns': df.columns.tolist(),
            'dtypes': df.dtypes.to_dict(),
            'memory_usage': df.memory_usage(deep=True).sum() / 1024,  # KB
        }
        
        print(f"📊 Catalog Info for '{catalog_name}':")
        print(f"   Rows: {info['num_rows']}")
        print(f"   Columns: {info['num_columns']}")
        print(f"   Column names: {', '.join(info['columns'])}")
        print(f"   Memory usage: {info['memory_usage']:.2f} KB")
        
        return info
    
    def save_filtered_catalog(self, catalog_name, output_path, **filters):
        """
        Save a filtered subset of a catalog to a new CSV file.
        
        Parameters
        ----------
        catalog_name : str
            Name of the source catalog file.
        output_path : str
            Path where the filtered catalog should be saved.
        **filters : keyword arguments
            Filtering conditions for the catalog.
        
        Returns
        -------
        pd.DataFrame
            The filtered DataFrame that was saved.
        """
        result = self.query_catalog(catalog_name, **filters)
        
        if result is None or len(result) == 0:
            print("⚠️  No rows matched the filter criteria. File not saved.")
            return None
        
        try:
            result.to_csv(output_path, index=False)
            print(f"✅ Filtered catalog saved to: {output_path}")
            return result
        except Exception as e:
            print(f"❌ Failed to save filtered catalog. Error: {e}")
            return None
    
    def display_catalog(self, catalog_name, max_rows=10):
        """
        Display a catalog in a formatted way.
        
        Parameters
        ----------
        catalog_name : str
            Name of the catalog file.
        max_rows : int, optional
            Maximum number of rows to display. Default is 10.
        """
        df = self.get_catalog(catalog_name)
        
        if df is None:
            return
        
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        pd.set_option('display.max_colwidth', None)
        
        print(f"\n📄 Displaying '{catalog_name}' (first {max_rows} rows):")
        print(df.head(max_rows))
        print(f"\nTotal rows: {len(df)}")
