


class ErrorTracker {

    constructor(errorDivID, successMsg) {
        this.hasError = false;
        this.message = null;
        this.successMsg = successMsg;
        this.target = errorDivID;
    }

    error(message) {
        this.hasError = true;
        this.message = message;

        const target = document.querySelector(`#${this.target}`);
        target.innerHTML = this.message;
        target.style.color = "red";
        target.style.display = "block";
    }

    clear() {
        this.hasError = false;
        this.message = this.successMsg;

        const target = document.querySelector(`#${this.target}`);
        target.innerHTML = this.message;
        target.style.color = "green";
        target.style.display = "block";
    }

}