document.addEventListener("DOMContentLoaded", function() {
    console.log("Website loaded successfully");
});
function toggleOlderNews() {
    var olderNews = document.getElementById("older-news");
    olderNews.style.display = olderNews.style.display === "none" ? "block" : "none";
}
document.addEventListener("DOMContentLoaded", function() {
    document.querySelectorAll("li").forEach(function(li) {
        li.innerHTML = li.innerHTML.replace(
            /(https?:\/\/[^\s<]+)/g,
            '<a href="$1" target="_blank" rel="noopener noreferrer">$1</a>'
        );
    });
});
