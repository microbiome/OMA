## Rendering the fennica project in quarto 

Execute the script `scriptConversionOMA.R` in the OMA directory.
To do it in terminal use `Rscript scriptConversionOMA.R`
The script will automatically convert the .Rmd files in quarto and then render the directory in the `_book` folder (can be changed in the `_quarto.yml` file).

To see what the book looks like, use the command below in the OMA director :

`quarto::quarto_preview()` in a R prompt 
or use :
`quarto preview` in the shell prompt


## Convert a Rmardown file to quarto format

In a R command prompt, use the following command : 
```R
knitr::convert_chunk_header("filename.Rmd", output = "filename.qmd")
``` 

The output is necessary to create a new file, otherwise the file converted in quarto is only printed. 


## Rendering a qmd document 

With the quarto package loaded, use in a shell prompt :
```sh
quarto render filename.qmd
```

or in a R prompt :
```R
quarto::quarto_render("filename.qmd")
```

-> The output of this command is set by default to create a html document of the same name.

## Renderin a quarto project

In the folder you want to render, use in a shell prompt :
```sh
quarto render 
```

-> all the qmd file in the current folder will be rendered. 

In the `_quarto.yml` file, you can precise what file you want to render and add treatment options.

---

Here are two links that can be helpfull : 
-> converting rmd to quarto for beginner
https://laderast.github.io/qmd_rmd/#/title-slide 

-> the official quarto documentation site
https://quarto.org/docs/guide/