result_data = "empty";

// To reset file input, not sure whether needed anymore
Shiny.addCustomMessageHandler("resetFileInputHandler", function(x) {      
        var id = "#" + x + "_progress";
        var idFile = "#" + x;
        var idBar = id + " .bar";
        $(id).css("visibility", "hidden");
        $(idBar).css("width", "0%");
        $(id).addClass("active");
        $(idFile).replaceWith(idFile = $(idFile).clone(true));
    });
    
// read data when sent to app    
$(document).on("shiny:connected", function() {
 window.addEventListener("message", displayMessage, false);
 function displayMessage(evt) { 
   var inmessage = JSON.parse(evt.data);
   console.log(inmessage); 
   console.log("read message");
   if (inmessage === "Retrieve results") {
     evt.source.postMessage("VSClust: result request received",evt.origin);
     setTimeout(evt.source.postMessage(result_data, evt.origin), 5000);
     Shiny.setInputValue("retrieve_output", "Get data");
   } else {
     evt.source.postMessage("VSClust: data received",evt.origin);
     Shiny.setInputValue("extdata", evt.data);
   }
 }
});


 // check ext window and retrieve results
shinyjs.send_results = function(params)
{
  
  result_data = params.dat;
}   
 
