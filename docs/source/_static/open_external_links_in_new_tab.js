document.addEventListener("DOMContentLoaded", function () {
    const links = document.querySelectorAll("a[href^='http']");
    for (const link of links) {
        link.setAttribute("target", "_blank");
    }
});