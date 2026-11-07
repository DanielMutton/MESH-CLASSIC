Configuring the model outputs via the CLASSIC code and Output Variable Editor (OVE) {#xmlSystem}
========

All CLASSIC model outputs are written to netcdf files. Outputs are first specified in an xml file where the relevant metadata is assigned. Then within the CLASSIC code the variable of interest is appropriately prepared for writing at the half-hourly, daily, monthly, or annual timestep (see @ref writeOutput).

The Output Variable Editor (OVE) is a web interface that allows a user to configure the output variables for CLASSIC. Using this interface the user can load an existing XML file and add/remove variables and variants. Once the changes are complete, the user may download an updated version of the XML configuration file.

# The XML Document Structure {#xmlStruct}
The XML document follows a hierarchical structure which has the **variableSet** as a root node. This node has a number of **group** child nodes, as well as attributes such as *type, version and created*. The **group** nodes have a number of **variable** child nodes.

The **variable** nodes feature several attributes (e.g. *id, includeBareGround*) and fields (e.g. *standardName, longName, shortName,* etc. They also have a collection of **variant** child nodes. These variants are the individual instances of the variable, which can be identified by the unique **nameInCode** identifier. They also have a specific **timeFrequency** and **outputForm**.

Example of xml:

        <?xml version="1.0"?>
        <variableSet type="CLASS" version="1.2" created="2018/8/7">
            <group type="class">
                <variable id="0" includeBareGround="false" dormant="false">
                    <standardName>Net Shortwave Surface Radiation</standardName>
                    <longName>Net downward shortwave radiation at the surface</longName>
                    <shortName>rss</shortName>
                    <defaultUnits>W/$m^2$</defaultUnits>
                    <variant>
                        <timeFrequency>annually</timeFrequency>
                        <outputForm>grid</outputForm>
                        <nameInCode>fsstar_yr</nameInCode>
                        <units>W/$m^2$</units>
                    </variant>
                    <variant>
                        <timeFrequency>monthly</timeFrequency>
                        <outputForm>grid</outputForm>
                        <nameInCode>fsstar_mo</nameInCode>
                        <units>W/$m^2$</units>
                    </variant>

                    ...

                  </variable>
              </group>
          </variableSet>

# Output Variable Editor (OVE)

## Opening Instructions
The *Output Variable Editor* folder is located in the *tools* folder of your local CLASS directory.
Open tools/outputVariableEditor/index.html in your favorite browser and edit the output variable descriptors using the web interface.

## Work flow
The various process steps are seen on the top of the page as tabs. The user first loads the existing XML document as a template, then optionally add/removes variables, then edits the variants associated with the variables, and finally generates the output XML that included the new variables and variant changes.
### Load the Existing XML
Browse to an existing XML file to add/remove variables and variants. Once the file has been selected, press the "Load XML" button to load the data. An existing file is provided in the configurationFiles folder. Choose the one with the newest version number.
### Add/Remove Variable
Fill in the **Standard Name**, **Short Name**, **Long Name**, **Units** fields and select a **Group** to classify the new variable. Select the **Include bare ground** check-box, if the variable includes bare ground. One may also use **Add variable group** to add a group or **Remove variable** to delete an existing variable.
### Variants editor
In the variants editor, the user can find a list of all of the existing variables. Click on one to find a detailed description of the variable as well as a list of potential variants. Check the box for a particular variant you're interested in and enter the desired **name in code**.

_**Warning**: Unchecking the box next to a variant and navigating to a different tab will remove any text entered in the **name in code** field_

You may also use the search bar at the top of the window to find a specific keyword in your list of variables, and/or use the **Make All Variables Dormant** button to make all currently-active variables inactive.

### Output XML
Once all the variable and variant changes are complete, proceed to the **Output XML** tab, where one can preview the newly generated XML document. Also, the website will trigger a download command, so the generated XML document will be pushed to the user for download.


# Validation / Schema
A relaxNG schema is provided for validation purposes.
In order to validate the XML Document, navigate to the schema folder and use the following command:

`java -jar jing.jar schema_1.0.xml path/myFile.xml`

where schema_1.0.xml is the current version of the schema and path/myFile.xml is the XML Document file of your choosing.

The schema will produce warnings and errors if there is any problem with the XML file, but it will **not** produce any output at all if everything is alright.

# Manual editing
There are certain functionalities that the present interface is lacking, like editing existing variable names, removing groups, editing groups and so on.
These tasks can be easily achieved by manually editing the XML input file.

**Do not manually edit the id fields of the variables, as that can cause serious problems with the variables**

# Editing the CLASSIC code to write the output variable {#writeOutput}

Then within the CLASSIC code (@ref prepareOutputs.f90) the variable of interest is appropriately prepared for writing at the half-hourly (@ref prepareoutputs::class_hh_w), daily (@ref prepareoutputs::ctem_daily_aw, @ref prepareoutputs::class_daily_aw), monthly (@ref prepareoutputs::ctem_monthly_aw, @ref prepareoutputs::class_monthly_aw) or annual (@ref prepareoutputs::ctem_annual_aw, @ref prepareoutputs::class_annual_aw) timestep.

All writes to the netcdf file in @ref prepareOutputs.f90 use some form of the following:

        call writeOutput1D(lonLocalIndex,latLocalIndex,'laimaxg_mo' ,timeStamp,'lai', [laimaxg_mo(i,m,:)])

Here there are two important considerations. First is the **key** ('laimaxg_mo' in this example). The key is used to associate this write to the appropriate netcdf output file, called the nameInCode in the xml file. The **shortName** (here it is 'lai') is the name of the variable in the output netcdf file. You will see the nameInCode and shortName in the xml snippet @ref xmlStruct "above". So your nameInCode and shortName within prepareOutputs.f90 and in your xml file must match for the proper netcdf to match to the proper model variable to be output.
